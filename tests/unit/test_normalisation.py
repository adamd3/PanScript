#!/usr/bin/env python3
"""Unit tests for count merging and normalisation on hand-built synthetic data.

The synthetic dataset is small enough that every expected value is computed
independently here (in numpy) from the formulas in docs/normalisation.md, and
compared with what the pipeline scripts in bin/ produce.

Run from the repo root (needs python3 + pandas, and Rscript with DESeq2,
optparse, matrixStats and tidyverse, e.g. inside the pipeline container):

    python3 -m unittest discover -s tests/unit -v
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"

# --------------------------------------------------------------------------
# Synthetic dataset
# --------------------------------------------------------------------------
# Three isolates (dna ids dA, dB, dC) -> three RNA samples (sA, sB, sC).
# Gene groups:
#   core1..core12 : single locus in every isolate, lengths differ by isolate
#   acc1          : accessory, absent from isolate C
#   acc2          : accessory, present only in isolate B
#   split1        : split into 2 loci in isolate A, single locus in B and C
#   split2        : split in isolate A, but NEITHER locus was quantified (both
#                   missing from A's FASTA) -> must be NA in A, not 0
#   split3        : split in isolate A, only one of its two loci quantified
#                   (partial; flagged by the missing-loci report, PROGRESS P7)
# A locus whose est_counts is None is listed by Panaroo but absent from the
# kallisto output (i.e. missing from the isolate FASTA).
SAMPLES = ["sA", "sB", "sC"]
DNA = {"sA": "dA", "sB": "dB", "sC": "dC"}
N_CORE = 12


def build_dataset():
    rng = np.random.default_rng(42)
    genes = {}  # gene -> {sample: [(locus, est_counts, eff_length), ...] or None}
    depth = {"sA": 1.0, "sB": 2.0, "sC": 0.5}  # library-size multipliers
    for i in range(1, N_CORE + 1):
        g = f"core{i}"
        base_len = rng.integers(400, 3000)
        expr = rng.integers(20, 2000)  # reads per kb at depth 1
        genes[g] = {}
        for s in SAMPLES:
            L = float(base_len + rng.integers(-60, 60))  # ortholog length varies
            c = float(round(expr * depth[s] * L / 1000 * rng.uniform(0.8, 1.25)))
            genes[g][s] = [(f"{DNA[s]}_{g}", c, L)]
    genes["acc1"] = {"sA": [("dA_acc1", 300.0, 900.0)],
                     "sB": [("dB_acc1", 520.0, 880.0)], "sC": None}
    genes["acc2"] = {"sA": None, "sB": [("dB_acc2", 75.0, 610.0)], "sC": None}
    genes["split1"] = {"sA": [("dA_split1a", 120.0, 700.0), ("dA_split1b", 80.0, 500.0)],
                       "sB": [("dB_split1", 410.0, 1190.0)],
                       "sC": [("dC_split1", 95.0, 1210.0)]}
    genes["split2"] = {"sA": [("dA_split2a", None, None), ("dA_split2b", None, None)],
                       "sB": [("dB_split2", 260.0, 830.0)],
                       "sC": [("dC_split2", 70.0, 845.0)]}
    genes["split3"] = {"sA": [("dA_split3a", 55.0, 400.0), ("dA_split3b", None, None)],
                       "sB": [("dB_split3", 330.0, 1020.0)],
                       "sC": [("dC_split3", 81.0, 1000.0)]}
    return genes


def write_inputs(d: Path, genes):
    # Panaroo-style gene_presence_absence (comma-separated)
    rows = []
    for g, per in genes.items():
        row = {"Gene": g, "Non-unique Gene name": "", "Annotation": "",
               "REF": f"ref_{g}"}
        for s in SAMPLES:
            loci = per[s]
            row[DNA[s]] = ";".join(l[0] for l in loci) if loci else ""
        rows.append(row)
    pd.DataFrame(rows).to_csv(d / "gpa.csv", index=False)

    meta = pd.DataFrame({"rna_sample_id": [f"R{s}" for s in SAMPLES],
                         "dna_sample_id": [DNA[s] for s in SAMPLES],
                         "sample_name": SAMPLES})
    meta.to_csv(d / "meta.tsv", sep="\t", index=False)

    for s in SAMPLES:
        kd = d / f"kallisto_{s}"
        kd.mkdir()
        recs = [(l[0], l[2] + 150, l[2], l[1], 0.0)
                for per in genes.values() if per[s] for l in per[s]
                if l[1] is not None]
        pd.DataFrame(recs, columns=["target_id", "length", "eff_length",
                                    "est_counts", "tpm"]).to_csv(
            kd / "abundance.tsv", sep="\t", index=False)


# --------------------------------------------------------------------------
# Independent reference implementation (current pipeline behaviour)
# --------------------------------------------------------------------------
def ref_merged(genes, field):
    """Merged est_counts (field=1) or eff_length (field=2).

    Split loci are SUMMED over the loci that were quantified (P5); NA if the
    gene is absent or none of its loci were quantified.
    """
    out = pd.DataFrame(index=list(genes), columns=SAMPLES, dtype=float)
    for g, per in genes.items():
        for s in SAMPLES:
            vals = [l[field] for l in (per[s] or []) if l[field] is not None]
            if vals:
                out.loc[g, s] = np.sum(vals)
    return out


def ref_size_factors(raw):
    """DESeq2 median-of-ratios on rounded, complete-case (core) rows."""
    k = np.round(raw.dropna().to_numpy())
    with np.errstate(divide="ignore"):
        logk = np.log(k)
    loggeo = logk.mean(axis=1)
    keep = np.isfinite(loggeo)
    ratios = logk[keep] - loggeo[keep, None]
    return pd.Series(np.exp(np.median(ratios, axis=0)), index=raw.columns)


def ref_outputs(genes):
    counts = ref_merged(genes, 1)
    lens = ref_merged(genes, 2)
    med = lens.median(axis=1, skipna=True)
    raw = counts.div(lens).mul(med, axis=0)
    sf = ref_size_factors(raw)
    norm = raw.div(sf, axis=1)
    # nTPM (P2): own-length rate / size factor, scaled by one global constant
    rate = counts.div(lens).div(sf, axis=1)
    K = 1e6 / rate.sum(axis=0, skipna=True).mean()
    return {"counts": counts, "lens": lens, "raw": raw, "sf": sf,
            "norm": norm, "ntpm": rate * K}  # all linear scale (P6)


# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------
def read_matrix(p):
    return pd.read_csv(p, sep="\t", index_col=0).astype(float)


def run_scripts(genes):
    """Write inputs for `genes` to a temp dir and run the bin/ scripts on them."""
    if shutil.which("Rscript") is None:
        raise unittest.SkipTest("Rscript not available")
    tmp = Path(tempfile.mkdtemp(prefix="sq_unit_"))
    write_inputs(tmp, genes)
    env = dict(os.environ, PATH=f"{BIN}{os.pathsep}{os.environ['PATH']}")

    def run(cmd):
        r = subprocess.run(cmd, cwd=tmp, env=env, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"{cmd[0]} failed:\n{r.stdout}\n{r.stderr}")

    for script, out in [("merge_kallisto_counts.py", "counts.tsv"),
                        ("merge_kallisto_lens.py", "lens.tsv")]:
        run([sys.executable, str(BIN / script), "--gene_presence_absence=gpa.csv",
             "--metadata_merged=meta.tsv", f"--outf={out}"])
    (tmp / "out").mkdir()
    run(["Rscript", str(BIN / "DESeq2_normalise_counts.R"), "-c", "counts.tsv",
         "-l", "lens.tsv", "-o", "out"])
    return tmp


class NormalisationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.genes = build_dataset()
        cls.tmp = run_scripts(cls.genes)
        cls.expected = ref_outputs(cls.genes)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def assertFrameClose(self, got, exp, rtol=1e-9):
        got = got.loc[exp.index, exp.columns]
        np.testing.assert_array_equal(got.isna().to_numpy(), exp.isna().to_numpy(),
                                      err_msg="NA pattern differs")
        np.testing.assert_allclose(got.to_numpy(float), exp.to_numpy(float),
                                   rtol=rtol, equal_nan=True)

    def test_merged_counts(self):
        self.assertFrameClose(read_matrix(self.tmp / "counts.tsv"), self.expected["counts"])

    def test_merged_lengths(self):
        self.assertFrameClose(read_matrix(self.tmp / "lens.tsv"), self.expected["lens"])

    def test_split_gene_summed(self):
        counts = read_matrix(self.tmp / "counts.tsv")
        lens = read_matrix(self.tmp / "lens.tsv")
        self.assertEqual(counts.loc["split1", "sA"], 120.0 + 80.0)
        self.assertEqual(lens.loc["split1", "sA"], 700.0 + 500.0)
        # per-isolate rate equals the pooled rate over its loci
        raw = read_matrix(self.tmp / "out" / "raw_counts.tsv")
        med = np.median([1200.0, 1190.0, 1210.0])
        self.assertAlmostEqual(raw.loc["split1", "sA"], 200.0 / 1200.0 * med)

    def test_unquantified_split_gene_is_na_not_zero(self):
        for f in ["counts.tsv", "lens.tsv"]:
            self.assertTrue(np.isnan(read_matrix(self.tmp / f).loc["split2", "sA"]), f)

    def test_partially_quantified_split_gene_uses_present_loci(self):
        self.assertEqual(read_matrix(self.tmp / "counts.tsv").loc["split3", "sA"], 55.0)
        self.assertEqual(read_matrix(self.tmp / "lens.tsv").loc["split3", "sA"], 400.0)

    def test_accessory_genes_are_na_not_zero(self):
        raw = read_matrix(self.tmp / "out" / "raw_counts.tsv")
        self.assertTrue(np.isnan(raw.loc["acc1", "sC"]))
        self.assertTrue(raw.loc[["acc2"], ["sA", "sC"]].isna().all(axis=None))

    def test_raw_counts_length_scaled(self):
        self.assertFrameClose(read_matrix(self.tmp / "out" / "raw_counts.tsv"),
                              self.expected["raw"])

    def test_norm_counts(self):
        self.assertFrameClose(read_matrix(self.tmp / "out" / "norm_counts.tsv"),
                              self.expected["norm"], rtol=1e-7)

    def test_ntpm(self):
        self.assertFrameClose(read_matrix(self.tmp / "out" / "ntpm.tsv"),
                              self.expected["ntpm"], rtol=1e-7)

    def test_ntpm_properties(self):
        ntpm = read_matrix(self.tmp / "out" / "ntpm.tsv")
        counts = read_matrix(self.tmp / "counts.tsv")
        lens = read_matrix(self.tmp / "lens.tsv")
        sf = self.expected["sf"]
        # per-sample totals average exactly 1e6 (but are not each 1e6)
        self.assertAlmostEqual(ntpm.sum(axis=0).mean(), 1e6, delta=1e-3)
        # within a sample, nTPM is proportional to counts / own length
        # (i.e. comparable between genes, like TPM)
        for s in SAMPLES:
            r = (ntpm[s] / (counts[s] / lens[s])).dropna()
            self.assertLess(r.std() / r.mean(), 1e-9, s)
        # median length cancels: nTPM == norm_counts / median length * const
        norm = read_matrix(self.tmp / "out" / "norm_counts.tsv")
        med = lens.median(axis=1, skipna=True)
        r = (ntpm / norm.div(med, axis=0)).stack()
        self.assertLess(r.std() / r.mean(), 1e-9)
        # across samples, a core gene's nTPM ratio is its rate ratio divided by
        # the size-factor ratio (library size from size factors, not totals)
        g = "core1"
        rate = counts.loc[g] / lens.loc[g]
        np.testing.assert_allclose((ntpm.loc[g] / ntpm.loc[g, "sA"]).to_numpy(),
                                   ((rate / sf) / (rate["sA"] / sf["sA"])).to_numpy(),
                                   rtol=1e-6)

    def test_no_rpkm_output(self):
        self.assertFalse((self.tmp / "out" / "rpkm_counts.tsv").exists())

    def test_size_factors_track_depth(self):
        # depth multipliers 1 : 2 : 0.5 should be recovered approximately
        sf = self.expected["sf"]
        ratio = (sf / sf["sA"]).to_numpy()
        np.testing.assert_allclose(ratio, [1.0, 2.0, 0.5], rtol=0.15)


class AccessoryGeneRobustnessTest(unittest.TestCase):
    """Why nTPM rather than TPM (docs/normalisation.md, step 5).

    Give isolate B one extra, very highly expressed accessory gene, with
    everything else identical. Per-sample-total TPM then deflates every core
    gene in B. nTPM comparisons between samples must not change.
    """

    @classmethod
    def setUpClass(cls):
        base = build_dataset()
        big = build_dataset()
        total_b = sum(l[1] for per in big.values() if per["sB"]
                      for l in per["sB"] if l[1] is not None)
        big["accBig"] = {"sA": None, "sB": [("dB_accBig", total_b, 800.0)], "sC": None}
        cls.dirs = {"base": run_scripts(base), "big": run_scripts(big)}

    @classmethod
    def tearDownClass(cls):
        for d in cls.dirs.values():
            shutil.rmtree(d, ignore_errors=True)

    def core(self, which, name):
        m = read_matrix(self.dirs[which] / name)
        return m.loc[[f"core{i}" for i in range(1, N_CORE + 1)]]

    def test_ntpm_between_sample_comparisons_unchanged(self):
        base, big = self.core("base", "out/ntpm.tsv"), self.core("big", "out/ntpm.tsv")
        # ratios between samples for every core gene are identical
        np.testing.assert_allclose(big.div(big["sA"], axis=0).to_numpy(),
                                   base.div(base["sA"], axis=0).to_numpy(), rtol=1e-9)

    def test_naive_tpm_would_be_distorted(self):
        # demonstrates the problem nTPM avoids: per-sample-total TPM of core
        # genes in B roughly halves when B carries the big accessory gene
        def tpm(which):
            c = read_matrix(self.dirs[which] / "counts.tsv")
            l = read_matrix(self.dirs[which] / "lens.tsv")
            rate = c / l
            return (rate / rate.sum(axis=0) * 1e6).loc[[f"core{i}" for i in range(1, N_CORE + 1)]]
        shift = (tpm("big")["sB"] / tpm("base")["sB"]).median()
        self.assertLess(shift, 0.8)


if __name__ == "__main__":
    unittest.main()

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
                for per in genes.values() if per[s] for l in per[s]]
        pd.DataFrame(recs, columns=["target_id", "length", "eff_length",
                                    "est_counts", "tpm"]).to_csv(
            kd / "abundance.tsv", sep="\t", index=False)


# --------------------------------------------------------------------------
# Independent reference implementation (current pipeline behaviour)
# --------------------------------------------------------------------------
def ref_merged(genes, field):
    """Merged est_counts (field=1) or eff_length (field=2); split loci averaged."""
    out = pd.DataFrame(index=list(genes), columns=SAMPLES, dtype=float)
    for g, per in genes.items():
        for s in SAMPLES:
            if per[s]:
                out.loc[g, s] = np.mean([l[field] for l in per[s]])
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
    lib = raw.sum(axis=0, skipna=True) / 1e6
    rpkm = raw.div(lib, axis=1).div(med / 1e3, axis=0)
    return {"counts": counts, "lens": lens, "raw": raw, "sf": sf,
            "norm": np.log2(norm + 1), "rpkm": np.log2(rpkm + 1)}


# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------
def read_matrix(p):
    return pd.read_csv(p, sep="\t", index_col=0).astype(float)


class NormalisationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if shutil.which("Rscript") is None:
            raise unittest.SkipTest("Rscript not available")
        cls.tmp = Path(tempfile.mkdtemp(prefix="sq_unit_"))
        cls.genes = build_dataset()
        write_inputs(cls.tmp, cls.genes)
        env = dict(os.environ, PATH=f"{BIN}{os.pathsep}{os.environ['PATH']}")

        def run(cmd):
            r = subprocess.run(cmd, cwd=cls.tmp, env=env, capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError(f"{cmd[0]} failed:\n{r.stdout}\n{r.stderr}")

        for script, out in [("merge_kallisto_counts.py", "counts.tsv"),
                            ("merge_kallisto_lens.py", "lens.tsv")]:
            run([sys.executable, str(BIN / script), "--gene_presence_absence=gpa.csv",
                 "--metadata_merged=meta.tsv", f"--outf={out}"])
        # gene-subset file is required by the current CLI but unused with -p FALSE
        pd.DataFrame({"gene": list(cls.genes)}).to_csv(cls.tmp / "subset.tsv",
                                                       sep="\t", index=False)
        (cls.tmp / "out").mkdir()
        run(["Rscript", str(BIN / "DESeq2_normalise_counts.R"), "-c", "counts.tsv",
             "-l", "lens.tsv", "-g", "subset.tsv", "-p", "FALSE", "-t", "TRUE",
             "-o", "out"])
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

    def test_rpkm_counts(self):
        self.assertFrameClose(read_matrix(self.tmp / "out" / "rpkm_counts.tsv"),
                              self.expected["rpkm"], rtol=1e-7)

    def test_size_factors_track_depth(self):
        # depth multipliers 1 : 2 : 0.5 should be recovered approximately
        sf = self.expected["sf"]
        ratio = (sf / sf["sA"]).to_numpy()
        np.testing.assert_allclose(ratio, [1.0, 2.0, 0.5], rtol=0.15)


if __name__ == "__main__":
    unittest.main()

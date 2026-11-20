#!/usr/bin/env python3
"""Check a test_pe run against the simulation truth (test_data/pe/truth.tsv.gz).

Usage: python3 tests/check_pe_truth.py <results_dir>

Checks that:
  - every sample's strandedness was inferred as reverse (2), as simulated
  - kallisto's est_counts per locus track the simulated fragment counts
  - the merged counts per gene group track the truth summed per group
Exit code 0 if all pass, 1 otherwise.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
TRUTH = REPO / "test_data" / "pe" / "truth.tsv.gz"
MIN_R = 0.98          # Pearson r on log1p scale, loci/groups with >= 5 fragments
MAX_MEDIAN_REL_ERR = 0.10


def main(results):
    results = Path(results)
    truth = pd.read_csv(TRUTH, sep="\t")
    ok = True

    def report(name, passed, detail):
        nonlocal ok
        ok &= passed
        print(f"{'PASS' if passed else 'FAIL'}  {name}: {detail}")

    strand = pd.read_csv(results / "pipeline_info" / "strandedness_check.tsv", sep="\t")
    report("strandedness", (strand["used_strandedness"] == 2).all(),
           dict(zip(strand["sample"], strand["used_strandedness"])))

    merged = pd.read_csv(next((results / "kallisto_quant").glob("kallisto_merged_counts_*.tsv")),
                         sep="\t", index_col=0)
    for sample, t in truth.groupby("sample_name"):
        ab = pd.read_csv(results / "kallisto_quant" / f"kallisto_{sample}" / "abundance.tsv",
                         sep="\t", index_col=0)
        tl = t.set_index("locus")["fragments"].reindex(ab.index, fill_value=0)
        keep = tl >= 5
        r = np.corrcoef(np.log1p(tl[keep]), np.log1p(ab["est_counts"][keep]))[0, 1]
        rel = (np.abs(ab["est_counts"][keep] - tl[keep]) / tl[keep]).median()
        report(f"{sample} loci", r >= MIN_R and rel <= MAX_MEDIAN_REL_ERR,
               f"r = {r:.4f}, median relative error = {rel:.3f} ({keep.sum()} loci)")

        tg = t.groupby("gene_group")["fragments"].sum()
        got = merged[sample].reindex(tg.index)
        keep = (tg >= 5) & got.notna()
        r = np.corrcoef(np.log1p(tg[keep]), np.log1p(got[keep]))[0, 1]
        report(f"{sample} gene groups", r >= MIN_R,
               f"r = {r:.4f} ({keep.sum()} groups); groups NA in output: "
               f"{int(got.isna().sum())}")

    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1]))

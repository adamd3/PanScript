#!/usr/bin/env python3
"""Compare pipeline output matrices against the stored baseline.

Usage:
    python3 tests/compare_to_baseline.py <new_results_dir> [--baseline tests/baseline]
        [--transform FILE=TRANSFORM ...] [--rtol 1e-6] [--atol 1e-6]

For each baseline matrix, the matching file in the new results is located by
its stem (dynamic suffixes such as `_<gpa>_<meta>` are ignored), then compared
cell by cell. NA positions must match exactly; numeric values must agree within
tolerance. Differences are reported per file with the most-changed rows.

TRANSFORM is applied to the *baseline* values before comparing, so that an
intended change can be checked exactly, e.g. after removing the log2(x+1)
transform from outputs:
    --transform norm_counts=exp2m1      (baseline 2^x - 1 should equal new)

Exit code 0 if every compared file matches, 1 otherwise.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Matrices compared, keyed by a stable stem; value = glob pattern (relative to
# results dir for new runs, or to the baseline dir).
MATRICES = {
    "raw_counts": "gene_counts/raw_counts.tsv",
    "norm_counts": "gene_counts/norm_counts.tsv",
    "rpkm_counts": "gene_counts/rpkm_counts.tsv",
    "kallisto_merged_counts": "kallisto_quant/kallisto_merged_counts_*.tsv",
    "kallisto_merged_lens": "kallisto_quant/kallisto_merged_lens_*.tsv",
}

TRANSFORMS = {
    "identity": lambda x: x,
    "exp2m1": lambda x: np.power(2.0, x) - 1.0,
}


def find_one(root: Path, pattern: str):
    hits = sorted(root.glob(pattern))
    if len(hits) > 1:
        sys.exit(f"ERROR: pattern {pattern} matched several files in {root}: {hits}")
    return hits[0] if hits else None


def load(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df.apply(pd.to_numeric, errors="coerce")


def compare(name, base, new, rtol, atol, max_rows=10):
    ok = True
    msgs = []
    rows_only_base = base.index.difference(new.index)
    rows_only_new = new.index.difference(base.index)
    cols_only_base = base.columns.difference(new.columns)
    cols_only_new = new.columns.difference(base.columns)
    if len(rows_only_base) or len(rows_only_new):
        ok = False
        msgs.append(
            f"  rows: {len(rows_only_base)} only in baseline, "
            f"{len(rows_only_new)} only in new "
            f"(e.g. {list(rows_only_base[:3])} / {list(rows_only_new[:3])})"
        )
    if len(cols_only_base) or len(cols_only_new):
        ok = False
        msgs.append(
            f"  columns: only in baseline {list(cols_only_base)}, "
            f"only in new {list(cols_only_new)}"
        )

    rows = base.index.intersection(new.index)
    cols = base.columns.intersection(new.columns)
    b = base.loc[rows, cols].to_numpy(dtype=float)
    n = new.loc[rows, cols].to_numpy(dtype=float)

    na_b, na_n = np.isnan(b), np.isnan(n)
    na_mismatch = na_b != na_n
    if na_mismatch.any():
        ok = False
        msgs.append(f"  NA pattern differs in {int(na_mismatch.sum())} cells")

    both = ~na_b & ~na_n
    close = np.ones_like(b, dtype=bool)
    close[both] = np.isclose(b[both], n[both], rtol=rtol, atol=atol)
    n_diff = int((~close).sum())
    if n_diff:
        ok = False
        absdiff = np.where(both, np.abs(b - n), 0.0)
        row_max = absdiff.max(axis=1)
        order = np.argsort(-row_max)[:max_rows]
        msgs.append(
            f"  {n_diff} of {int(both.sum())} values differ "
            f"(max |diff| = {absdiff.max():.6g}); "
            f"{int((~close).any(axis=1).sum())} rows affected. Top rows:"
        )
        for i in order:
            if row_max[i] == 0:
                break
            msgs.append(f"    {rows[i]}: max |diff| {row_max[i]:.6g}")

    status = "MATCH" if ok else "DIFFERS"
    print(f"[{status}] {name}: {len(rows)} rows x {len(cols)} cols compared")
    for m in msgs:
        print(m)
    return ok


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("new_results", type=Path)
    p.add_argument("--baseline", type=Path,
                   default=Path(__file__).resolve().parent / "baseline")
    p.add_argument("--transform", action="append", default=[],
                   help="FILE=TRANSFORM applied to baseline values "
                        f"(transforms: {', '.join(TRANSFORMS)})")
    p.add_argument("--only", nargs="*", help="compare only these matrices")
    p.add_argument("--rtol", type=float, default=1e-6)
    p.add_argument("--atol", type=float, default=1e-6)
    args = p.parse_args()

    transforms = {}
    for t in args.transform:
        key, _, fn = t.partition("=")
        if key not in MATRICES or fn not in TRANSFORMS:
            sys.exit(f"ERROR: bad --transform {t}")
        transforms[key] = TRANSFORMS[fn]

    all_ok = True
    for key, pattern in MATRICES.items():
        if args.only and key not in args.only:
            continue
        bpath = find_one(args.baseline, pattern)
        npath = find_one(args.new_results, pattern)
        if bpath is None:
            print(f"[SKIP] {key}: not in baseline")
            continue
        if npath is None:
            print(f"[MISSING] {key}: not found in new results ({pattern})")
            all_ok = False
            continue
        base = load(bpath)
        base = base.apply(transforms.get(key, TRANSFORMS["identity"]))
        all_ok &= compare(key, base, load(npath), args.rtol, args.atol)

    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Check library strandedness against the --strandedness setting.

Two subcommands, used by the CHECK_STRANDEDNESS process:

  subsample  write the first N reads of one or two FASTQ files (gzip or plain)
  infer      compare kallisto pseudoalignment of the subsample run with
             --fr-stranded and with --rf-stranded, infer the library type and
             report whether it agrees with the --strandedness parameter

Inference (thresholds follow nf-core/rnaseq): with f = fr / (fr + rf),
  f >= 0.8         -> forward-stranded  (--strandedness 1)
  f <= 0.2         -> reverse-stranded  (--strandedness 2)
  0.4 <= f <= 0.6  -> unstranded        (--strandedness 0)
  otherwise        -> undetermined
"""

import argparse
import gzip
import json
import sys

LABELS = {0: "unstranded", 1: "forward", 2: "reverse"}


def open_any(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def subsample(args):
    for i, path in enumerate(args.reads, start=1):
        n = 0
        with open_any(path) as fin, gzip.open(f"{args.out}_{i}.fq.gz", "wt") as fout:
            for line_no, line in enumerate(fin):
                if line_no >= 4 * args.n:
                    break
                fout.write(line)
                n = line_no // 4 + 1
        print(f"{path}: wrote {n} reads", file=sys.stderr)


def classify(fr, rf, stranded=0.8, unstranded_band=(0.4, 0.6)):
    total = fr + rf
    if total == 0:
        return None, float("nan")
    f = fr / total
    if f >= stranded:
        return 1, f
    if rf / total >= stranded:   # not f <= 1 - stranded: 1 - 0.8 != 0.2 in floats
        return 2, f
    if unstranded_band[0] <= f <= unstranded_band[1]:
        return 0, f
    return None, f


def infer(args):
    fr = json.load(open(args.fr))
    rf = json.load(open(args.rf))
    inferred, f = classify(fr["n_pseudoaligned"], rf["n_pseudoaligned"])
    specified = args.specified
    agrees = inferred == specified if inferred is not None else None
    row = {
        "sample": args.sample,
        "reads_tested": fr["n_processed"],
        "pseudoaligned_fr": fr["n_pseudoaligned"],
        "pseudoaligned_rf": rf["n_pseudoaligned"],
        "fraction_fr": round(f, 4),
        "inferred": LABELS.get(inferred, "undetermined"),
        "inferred_strandedness": "" if inferred is None else inferred,
        "specified_strandedness": specified,
        "agrees": {True: "yes", False: "NO", None: "unknown"}[agrees],
    }
    with open(args.out, "w") as fh:
        fh.write("\t".join(row) + "\n")
        fh.write("\t".join(str(v) for v in row.values()) + "\n")

    pct = f"{f:.0%} of stranded pseudoalignments are forward"
    if agrees is False:
        print(f"Sample {args.sample}: library looks {LABELS[inferred]}-stranded ({pct}), but "
              f"--strandedness {specified} ({LABELS[specified]}) was used. Consider "
              f"--strandedness {inferred}. See pipeline_info/strandedness_check.tsv")
    elif agrees is None:
        print(f"Sample {args.sample}: could not determine library strandedness ({pct}). "
              f"Check --strandedness {specified} ({LABELS[specified]}). "
              f"See pipeline_info/strandedness_check.tsv")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("subsample")
    s.add_argument("--n", type=int, default=200000)
    s.add_argument("--out", required=True, help="prefix; writes <out>_1.fq.gz [<out>_2.fq.gz]")
    s.add_argument("reads", nargs="+")
    i = sub.add_parser("infer")
    i.add_argument("--fr", required=True, help="run_info.json from --fr-stranded run")
    i.add_argument("--rf", required=True, help="run_info.json from --rf-stranded run")
    i.add_argument("--sample", required=True)
    i.add_argument("--specified", type=int, choices=[0, 1, 2], required=True)
    i.add_argument("--out", required=True)
    a = p.parse_args()
    {"subsample": subsample, "infer": infer}[a.cmd](a)


if __name__ == "__main__":
    main()

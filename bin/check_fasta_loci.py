#!/usr/bin/env python3
"""Check that the loci Panaroo lists for an isolate are present in its FASTA.

Loci listed in the gene_presence_absence table but missing from the isolate's
gene FASTA cannot be quantified by kallisto, so their gene groups end up NA
(treated as absent) for that sample. This script reports them so the user knows
which genes were not quantified, and stops the run if most loci are missing,
which almost always means the wrong FASTA or a mismatched dna_sample_id.

Panaroo's status suffixes on locus IDs (_len, _pseudo, _stop) are not part of
the sequence ID, so a locus counts as present if its raw ID, or its ID with
those suffixes stripped, is in the FASTA (the merge scripts match the same way).

Outputs:
  --outf     one row per missing locus (sample, gene group, locus, flags)
  --summary  one row per sample with counts
  stdout     a one-line warning if any loci are missing (shown by Nextflow)

Flags:
  refound_stop       locus is a Panaroo "refound" gene with a premature stop
                     codon (*refound*_stop), typically a pseudogene that
                     annotation tools do not output as a CDS
  refound_len        Panaroo "refound" gene flagged with a length problem
                     (*refound*_len, e.g. length not a multiple of three)
  split_partial      the gene group is split into several loci in this isolate
                     and only some of them are missing; its counts and length
                     are summed over the loci present, so it is undercounted
  split_all_missing  split gene group with every locus missing (NA)
"""

import argparse
import gzip
import re
import sys

import pandas as pd

REFOUND_FLAGGED = re.compile(r"refound.*_(stop|len)$")

# Panaroo appends status suffixes to locus IDs in gene_presence_absence.csv
# (_len, _pseudo, _stop; sometimes several, e.g. _len_stop). They are not part
# of the sequence ID, so a locus is matched on its raw ID first and then with
# suffixes stripped one at a time, exactly as in merge_kallisto_*.py
# (PROGRESS.md R12).
PANAROO_SUFFIX = re.compile(r"_(len|pseudo|stop)$")


def resolve_locus(locus, ids):
    """Return the ID in `ids` that `locus` refers to, or `locus` if none does."""
    cur = locus
    while cur not in ids:
        stripped = PANAROO_SUFFIX.sub("", cur)
        if stripped == cur:
            return locus
        cur = stripped
    return cur


def fasta_ids(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    ids = set()
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gene_presence_absence", required=True)
    p.add_argument("--dna_sample_id", required=True)
    p.add_argument("--sample_name", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--max_missing_frac", type=float, default=0.5)
    p.add_argument("--outf", required=True)
    p.add_argument("--summary", required=True)
    a = p.parse_args()

    # read only the columns needed: the table can have thousands of isolates
    header = pd.read_csv(a.gene_presence_absence, nrows=0).columns
    if a.dna_sample_id not in header:
        sys.exit(f"ERROR: dna_sample_id '{a.dna_sample_id}' (sample {a.sample_name}) is not a "
                 f"column of the gene presence/absence file.")
    gpa = pd.read_csv(a.gene_presence_absence, low_memory=False, dtype=str,
                      usecols=["Gene", a.dna_sample_id])
    cells = gpa.set_index("Gene")[a.dna_sample_id].dropna()
    cells = cells[cells.str.strip() != ""]
    ids = fasta_ids(a.fasta)

    rows = []
    n_loci = 0
    groups_na = 0
    for gene, cell in cells.items():
        loci = [l.strip() for l in cell.split(";") if l.strip()]
        n_loci += len(loci)
        missing = [l for l in loci if resolve_locus(l, ids) not in ids]
        if not missing:
            continue
        split = len(loci) > 1
        all_missing = len(missing) == len(loci)
        groups_na += all_missing
        for locus in missing:
            flags = []
            m = REFOUND_FLAGGED.search(locus)
            if m:
                flags.append(f"refound_{m.group(1)}")
            if split:
                flags.append("split_all_missing" if all_missing else "split_partial")
            rows.append({"sample": a.sample_name, "dna_sample_id": a.dna_sample_id,
                         "gene_group": gene, "locus": locus, "n_loci_in_group": len(loci),
                         "flags": ",".join(flags)})

    cols = ["sample", "dna_sample_id", "gene_group", "locus", "n_loci_in_group", "flags"]
    miss = pd.DataFrame(rows, columns=cols)
    miss.to_csv(a.outf, sep="\t", index=False)

    n_missing = len(miss)
    frac = n_missing / n_loci if n_loci else 0.0
    summary = {
        "sample": a.sample_name,
        "dna_sample_id": a.dna_sample_id,
        "gene_groups_listed": len(cells),
        "loci_listed": n_loci,
        "loci_missing_from_fasta": n_missing,
        "frac_loci_missing": round(frac, 6),
        "gene_groups_affected": miss["gene_group"].nunique(),
        "gene_groups_na_due_to_missing": groups_na,
        "loci_refound_stop": int(miss["flags"].str.contains("refound_stop").sum()),
        "loci_refound_len": int(miss["flags"].str.contains("refound_len").sum()),
        "gene_groups_split_partial": miss.loc[miss["flags"].str.contains("split_partial"),
                                              "gene_group"].nunique(),
    }
    pd.DataFrame([summary]).to_csv(a.summary, sep="\t", index=False)

    if n_loci == 0:
        sys.exit(f"ERROR: no loci listed for '{a.dna_sample_id}' (sample {a.sample_name}) in the "
                 f"gene presence/absence file.")
    if frac > a.max_missing_frac:
        sys.exit(f"ERROR: sample {a.sample_name}: {n_missing} of {n_loci} loci ({frac:.1%}) listed "
                 f"for '{a.dna_sample_id}' in the gene presence/absence file are missing from its "
                 f"FASTA ({a.fasta}), above the limit of {a.max_missing_frac:.1%} "
                 f"(--max_missing_loci). Check that the FASTA and dna_sample_id belong to the "
                 f"same isolate.")
    if n_missing:
        print(f"Sample {a.sample_name}: {n_missing} of {n_loci} loci ({frac:.1%}) listed by Panaroo "
              f"are missing from the FASTA ({summary['loci_refound_stop']} refound_stop, "
              f"{summary['loci_refound_len']} refound_len, "
              f"{summary['gene_groups_split_partial']} partially missing split genes); "
              f"{summary['gene_groups_na_due_to_missing']} gene groups will be NA for this sample. "
              f"See pipeline_info/missing_loci.tsv")


if __name__ == "__main__":
    main()

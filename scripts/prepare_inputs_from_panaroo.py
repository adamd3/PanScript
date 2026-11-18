#!/usr/bin/env python3
"""Prepare pipeline inputs from a Panaroo output directory.

Writes, for every isolate (genome) in the Panaroo run:
  <outdir>/fasta/<isolate>.fna   the isolate's gene (CDS) nucleotide sequences,
                                 headed by the same locus IDs Panaroo uses in
                                 gene_presence_absence.csv
and a metadata template:
  <outdir>/metadata_template.tsv one row per isolate with dna_sample_id and
                                 fasta filled in; add rna_sample_id,
                                 sample_name, fastq1 (and fastq2), and repeat a
                                 row for each RNA-seq sample of that isolate

Uses Panaroo's gene_data.csv, which contains every annotated gene and every
"refound" gene with its DNA sequence, so genes that annotation tools do not
output as CDS (for example refound pseudogenes) are included.

Panaroo appends status suffixes to locus IDs in gene_presence_absence.csv
(e.g. "_len", "_pseudo", or "_stop" in older versions) that are not part of
the sequence ID; these are removed when matching.

Usage:
  prepare_inputs_from_panaroo.py --panaroo_dir <panaroo_output> --outdir <dir>
      [--isolates ID1 ID2 ...]
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

STATUS_SUFFIX = re.compile(r"(_(?:len|pseudo|stop))+$")


def strip_status(locus):
    return STATUS_SUFFIX.sub("", locus)


def read_gene_data(path):
    """Yield (isolate, annotation_id, dna_sequence) from Panaroo gene_data.csv.

    Columns: gff_file, scaffold_name, clustering_id, annotation_id,
    prot_sequence, dna_sequence, gene_name, description. Panaroo writes rows
    by joining fields with commas without quoting, and the last field
    (description) can itself contain commas, so split at most 7 times.
    """
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split(",")
        expected = ["gff_file", "scaffold_name", "clustering_id", "annotation_id",
                    "prot_sequence", "dna_sequence", "gene_name", "description"]
        if header[:8] != expected:
            sys.exit(f"ERROR: unexpected gene_data.csv header: {header}")
        for line in fh:
            f = line.rstrip("\n").split(",", 7)
            if len(f) < 6:
                continue
            yield f[0], f[3], f[5]


def read_gpa_loci(path):
    """Return (isolates, {isolate: set of locus IDs, suffixes stripped})."""
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        # isolate columns: after the annotation columns. Panaroo's
        # gene_presence_absence.csv has 3 (Gene, Non-unique Gene name,
        # Annotation); the Roary-style file has 14.
        n_meta = 14 if "No. isolates" in header else 3
        isolates = header[n_meta:]
        loci = defaultdict(set)
        for row in reader:
            for iso, cell in zip(isolates, row[n_meta:]):
                for l in cell.split(";"):
                    if l.strip():
                        loci[iso].add(strip_status(l.strip()))
    return isolates, loci


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--panaroo_dir", required=True,
                   help="Panaroo output directory (needs gene_data.csv and "
                        "gene_presence_absence.csv)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--isolates", nargs="*",
                   help="only these isolates (default: all in the Panaroo run)")
    a = p.parse_args()

    gene_data = os.path.join(a.panaroo_dir, "gene_data.csv")
    gpa = os.path.join(a.panaroo_dir, "gene_presence_absence.csv")
    for f in (gene_data, gpa):
        if not os.path.exists(f):
            sys.exit(f"ERROR: {f} not found")

    isolates, gpa_loci = read_gpa_loci(gpa)
    if a.isolates:
        unknown = set(a.isolates) - set(isolates)
        if unknown:
            sys.exit(f"ERROR: not isolates in {gpa}: {sorted(unknown)}")
        isolates = [i for i in isolates if i in set(a.isolates)]

    seqs = defaultdict(dict)
    wanted = set(isolates)
    for iso, locus, dna in read_gene_data(gene_data):
        if iso in wanted and dna:
            seqs[iso][locus] = dna

    fasta_dir = os.path.join(a.outdir, "fasta")
    os.makedirs(fasta_dir, exist_ok=True)
    rows = []
    problems = 0
    for iso in isolates:
        if not seqs[iso]:
            print(f"WARNING: no sequences for isolate '{iso}' in gene_data.csv", file=sys.stderr)
            problems += 1
            continue
        out = os.path.join(fasta_dir, f"{iso}.fna")
        with open(out, "w") as fh:
            for locus, dna in seqs[iso].items():
                fh.write(f">{locus}\n")
                for i in range(0, len(dna), 60):
                    fh.write(dna[i:i + 60] + "\n")
        missing = gpa_loci[iso] - set(seqs[iso])
        msg = (f"{iso}: {len(seqs[iso])} sequences written; "
               f"{len(gpa_loci[iso]) - len(missing)}/{len(gpa_loci[iso])} loci in "
               f"gene_presence_absence.csv have a sequence")
        if missing:
            problems += 1
            msg += f" (missing, e.g.: {sorted(missing)[:3]})"
        print(msg, file=sys.stderr)
        rows.append({"rna_sample_id": "", "dna_sample_id": iso, "sample_name": "",
                     "fastq1": "", "fastq2": "", "fasta": os.path.abspath(out)})

    with open(os.path.join(a.outdir, "metadata_template.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["rna_sample_id", "dna_sample_id", "sample_name",
                                           "fastq1", "fastq2", "fasta"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} FASTA files to {fasta_dir} and "
          f"{os.path.join(a.outdir, 'metadata_template.tsv')}", file=sys.stderr)
    if problems:
        print(f"WARNING: {problems} isolate(s) had missing sequences; see messages above.",
              file=sys.stderr)


if __name__ == "__main__":
    main()

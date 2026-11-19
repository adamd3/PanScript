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

Either file may be gzip-compressed (gene_data.csv.gz,
gene_presence_absence.csv.gz). gene_data.csv can be very large (tens of GB for
hundreds of genomes); it is streamed line by line, so memory use depends only
on the sequences of the isolates selected (use --isolates to limit them).

Usage:
  prepare_inputs_from_panaroo.py --panaroo_dir <panaroo_output> --outdir <dir>
      [--isolates ID1 ID2 ...]
"""

import argparse
import csv
import gzip
import os
import re
import sys
from collections import defaultdict

STATUS_SUFFIX = re.compile(r"(_(?:len|pseudo|stop))+$")


def strip_status(locus):
    return STATUS_SUFFIX.sub("", locus)


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, newline="")


def find_input(panaroo_dir, name):
    """Return the path of <name> or <name>.gz in panaroo_dir (exits if neither)."""
    for candidate in (name, name + ".gz"):
        path = os.path.join(panaroo_dir, candidate)
        if os.path.exists(path):
            return path
    sys.exit(f"ERROR: {os.path.join(panaroo_dir, name)}(.gz) not found")


def read_gene_data(path, wanted=None):
    """Yield (isolate, annotation_id, dna_sequence) from Panaroo gene_data.csv.

    Columns: gff_file, scaffold_name, clustering_id, annotation_id,
    prot_sequence, dna_sequence, gene_name, description. Panaroo writes rows
    by joining fields with commas without quoting, and the last field
    (description) can itself contain commas, so split at most 7 times.
    Rows for isolates not in `wanted` (if given) are skipped cheaply.
    """
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split(",")
        expected = ["gff_file", "scaffold_name", "clustering_id", "annotation_id",
                    "prot_sequence", "dna_sequence", "gene_name", "description"]
        if header[:8] != expected:
            sys.exit(f"ERROR: unexpected gene_data.csv header: {header}")
        for line in fh:
            if wanted is not None and line[:line.find(",")] not in wanted:
                continue
            f = line.rstrip("\n").split(",", 7)
            if len(f) < 6:
                continue
            yield f[0], f[3], f[5]


def read_gpa_isolates(path):
    """Return the isolate columns of gene_presence_absence.csv and how many
    annotation columns precede them (Panaroo's file has 3: Gene, Non-unique
    Gene name, Annotation; the Roary-style file has 14)."""
    with open_text(path) as fh:
        header = next(csv.reader(fh))
    n_meta = 14 if "No. isolates" in header else 3
    return header[n_meta:], n_meta


def read_gpa_loci(path, wanted):
    """Return {isolate: set of locus IDs, suffixes stripped} for `wanted`.

    Only the selected isolates are kept, so memory doesn't grow with the
    number of genomes in the Panaroo run."""
    isolates, n_meta = read_gpa_isolates(path)
    keep = [(i, iso) for i, iso in enumerate(isolates, start=n_meta) if iso in wanted]
    with open_text(path) as fh:
        reader = csv.reader(fh)
        next(reader)
        loci = defaultdict(set)
        for row in reader:
            for i, iso in keep:
                for l in row[i].split(";"):
                    if l.strip():
                        loci[iso].add(strip_status(l.strip()))
    return loci


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--panaroo_dir", required=True,
                   help="Panaroo output directory (needs gene_data.csv and "
                        "gene_presence_absence.csv, either optionally gzipped)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--isolates", nargs="*",
                   help="only these isolates (default: all in the Panaroo run)")
    a = p.parse_args()

    gene_data = find_input(a.panaroo_dir, "gene_data.csv")
    gpa = find_input(a.panaroo_dir, "gene_presence_absence.csv")

    isolates, _ = read_gpa_isolates(gpa)
    if a.isolates:
        unknown = set(a.isolates) - set(isolates)
        if unknown:
            sys.exit(f"ERROR: not isolates in {gpa}: {sorted(unknown)}")
        isolates = [i for i in isolates if i in set(a.isolates)]

    wanted = set(isolates)
    gpa_loci = read_gpa_loci(gpa, wanted)
    seqs = defaultdict(dict)
    for iso, locus, dna in read_gene_data(gene_data, wanted):
        if dna:
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

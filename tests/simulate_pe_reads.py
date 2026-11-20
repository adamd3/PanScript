#!/usr/bin/env python3
"""Simulate the paired-end test dataset (test_data/pe/), with known truth.

Reads are simulated from the three test isolates' own gene FASTAs
(test_data/SRR87372*.fna), so the PE test uses the same Panaroo table as the
single-end test. Deterministic (fixed seed, gzip mtime 0): re-running this
script reproduces the committed files byte for byte.

Library model (Illumina, dUTP-style reverse-stranded, as most bacterial kits):
  - each gene group has an expression level shared by all isolates
    (log-normal), times a per-isolate factor; a few groups are 8x up in one
    isolate; the replicate sample shares its isolate's levels plus small noise
  - fragments are drawn in proportion to expression x length; fragment length
    ~ Normal(200, 50), clipped to [60, 500] and to the gene length
  - 2 x 100 bp reads; read 1 is antisense (reverse complement of the
    fragment's 3' end), read 2 is sense, i.e. kallisto --rf-stranded
  - fragments shorter than 100 bp read into the TruSeq adapters (exercises
    fastp's paired-end adapter trimming)
  - substitution errors at 0.2% per base, a little higher towards read ends
  - samples differ in depth (exercises the size factors)

Writes test_data/pe/: <sample>_R{1,2}.fq.gz, meta_pe.tsv, truth.tsv.gz
(fragments simulated per locus and sample).

Usage (from the repo root): python3 tests/simulate_pe_reads.py
"""

import gzip
import io
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "test_data"
OUT = DATA / "pe"
SEED = 20261120
READ_LEN = 100
FRAG_MEAN, FRAG_SD, FRAG_MIN, FRAG_MAX = 200, 50, 60, 500
ERR = 0.002
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  # TruSeq, read 1 read-through
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"   # TruSeq, read 2 read-through

# sample_name, isolate (dna_sample_id), read pairs
SAMPLES = [("PE_A", "SRR8737281", 40000),
           ("PE_B", "SRR8737282", 60000),
           ("PE_C1", "SRR8737283", 30000),
           ("PE_C2", "SRR8737283", 45000)]   # replicate of isolate C

COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def locus_groups(isolate):
    """locus ID -> gene group, from the Panaroo table (suffixes stripped)."""
    gpa = pd.read_csv(DATA / "gene_presence_absence.tsv", usecols=["Gene", isolate],
                      dtype=str).dropna()
    out = {}
    for g, cell in zip(gpa["Gene"], gpa[isolate]):
        for l in cell.split(";"):
            for sfx in ("_len", "_pseudo", "_stop"):
                l = l[:-len(sfx)] if l.endswith(sfx) else l
            out[l] = g
    return out


def mutate(seq, rng):
    pos = np.arange(len(seq))
    p = ERR * (1 + 2 * pos / max(len(seq), 1))   # more errors towards the 3' end
    hit = rng.random(len(seq)) < p
    if not hit.any():
        return seq
    s = list(seq)
    for i in np.flatnonzero(hit):
        s[i] = "ACGT"[(("ACGT".find(s[i]) + rng.integers(1, 4)) % 4)]
    return "".join(s)


def quals(n):
    return "F" * max(n - 15, 0) + ":" * min(n, 15)


def main():
    rng = np.random.default_rng(SEED)
    OUT.mkdir(exist_ok=True)

    isolates = sorted({iso for _, iso, _ in SAMPLES})
    seqs = {iso: read_fasta(DATA / f"{iso}.fna") for iso in isolates}
    groups = {iso: locus_groups(iso) for iso in isolates}

    # expression per gene group (shared), per-isolate factors, a few DE groups
    all_groups = sorted({g for iso in isolates for g in groups[iso].values()})
    base = dict(zip(all_groups, rng.lognormal(0, 1.8, len(all_groups))))
    iso_level = {}
    for iso in isolates:
        lv = {}
        de = set(rng.choice(all_groups, 50, replace=False))
        for l, g in groups[iso].items():
            if l in seqs[iso] and len(seqs[iso][l]) >= FRAG_MIN:
                lv[l] = base[g] * rng.lognormal(0, 0.3) * (8 if g in de else 1)
        iso_level[iso] = lv

    meta_rows, truth_rows = [], []
    for name, iso, n_pairs in SAMPLES:
        loci = sorted(iso_level[iso])
        lens = np.array([len(seqs[iso][l]) for l in loci], float)
        lv = np.array([iso_level[iso][l] for l in loci]) * rng.lognormal(0, 0.1, len(loci))
        w = lv * lens
        n_frag = rng.multinomial(n_pairs, w / w.sum())
        r1 = io.StringIO()
        r2 = io.StringIO()
        k = 0
        for l, n in zip(loci, n_frag):
            if n == 0:
                continue
            truth_rows.append((name, iso, l, groups[iso][l], int(n)))
            seq = seqs[iso][l]
            for _ in range(n):
                flen = int(np.clip(round(rng.normal(FRAG_MEAN, FRAG_SD)), FRAG_MIN,
                                   min(FRAG_MAX, len(seq))))
                start = int(rng.integers(0, len(seq) - flen + 1))
                frag = seq[start:start + flen]           # sense strand of the RNA
                # read 1: antisense (from the fragment's 3' end); read 2: sense
                s1 = (revcomp(frag) + ADAPTER_R1 + "A" * READ_LEN)[:READ_LEN]
                s2 = (frag + ADAPTER_R2 + "A" * READ_LEN)[:READ_LEN]
                s1, s2 = mutate(s1, rng), mutate(s2, rng)
                k += 1
                rid = f"@{name}.{k}"
                r1.write(f"{rid}/1\n{s1}\n+\n{quals(READ_LEN)}\n")
                r2.write(f"{rid}/2\n{s2}\n+\n{quals(READ_LEN)}\n")
        for mate, buf in (("R1", r1), ("R2", r2)):
            # mtime=0 so the gzip bytes are reproducible
            with open(OUT / f"{name}_{mate}.fq.gz", "wb") as fh, \
                    gzip.GzipFile(fileobj=fh, mode="wb", mtime=0, filename="") as gz:
                gz.write(buf.getvalue().encode())
        meta_rows.append({"rna_sample_id": f"sim_{name}", "dna_sample_id": iso,
                          "sample_name": name, "isolate": iso,
                          "fastq1": f"{name}_R1.fq.gz", "fastq2": f"{name}_R2.fq.gz",
                          "fasta": f"../{iso}.fna"})
        print(f"{name}: {k} read pairs from {len(loci)} loci of {iso}")

    pd.DataFrame(meta_rows).to_csv(OUT / "meta_pe.tsv", sep="\t", index=False)
    pd.DataFrame(truth_rows, columns=["sample_name", "dna_sample_id", "locus", "gene_group",
                                      "fragments"]).to_csv(OUT / "truth.tsv.gz", sep="\t",
                                                           index=False, compression={"method": "gzip", "mtime": 0})


if __name__ == "__main__":
    main()

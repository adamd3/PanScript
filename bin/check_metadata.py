#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import os.path
import sys

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "sample_file",
        help = "tsv file containing sample metadata"
    )
    parser.add_argument("outf", help="File for output")
    args = parser.parse_args()
    check_meta(**vars(args))

def check_meta(sample_file, outf):
    sample_dat = pd.read_csv(sample_file, sep = "\t")
    # remove empty rows
    sample_dat.replace("", float("NaN"), inplace=True)
    sample_dat.dropna(subset = ["dna_sample_id"], inplace=True)
    # single-end and paired-end samples can't be mixed in one run: kallisto is
    # run in one mode for every sample (PROGRESS.md R5/Q4)
    has_fq2 = sample_dat["fastq2"].astype(str).str.strip().ne("") & sample_dat["fastq2"].notna()
    if has_fq2.any() and not has_fq2.all():
        single = sample_dat.loc[~has_fq2, "sample_name"].astype(str).tolist()
        paired = sample_dat.loc[has_fq2, "sample_name"].astype(str).tolist()
        sys.exit(
            "ERROR: the metadata mixes single-end samples (no fastq2: "
            + ", ".join(single) + ") and paired-end samples (fastq2 given: "
            + ", ".join(paired) + "). All samples in a run must be either single-end or "
            "paired-end. Run them as separate pipeline runs.")
    # sample names label every output (kallisto directories, matrix columns)
    dup = sample_dat["sample_name"][sample_dat["sample_name"].duplicated()].astype(str).unique()
    if len(dup):
        sys.exit("ERROR: sample_name must be unique; repeated: " + ", ".join(dup))
    # samples from the same isolate (dna_sample_id) share one kallisto index,
    # so they must use the same FASTA
    n_fasta = sample_dat.groupby("dna_sample_id")["fasta"].nunique()
    if (n_fasta > 1).any():
        bad = n_fasta[n_fasta > 1].index.astype(str).tolist()
        sys.exit("ERROR: each dna_sample_id must have a single fasta, but these have "
                 "more than one: " + ", ".join(bad) + ". Samples from the same isolate "
                 "must point to the same FASTA file.")
    sample_dat['paired'] = "1" if has_fq2.all() else "0"
    sample_dat.to_csv(outf, index=False, sep='\t')

if __name__ == "__main__":
    parse()

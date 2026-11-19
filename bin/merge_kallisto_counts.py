#!/usr/bin/env python3

import argparse
import pandas as pd
import os.path
import re


# Panaroo appends status suffixes to locus IDs in gene_presence_absence.csv
# (_len: length deviates from the group; _pseudo: invalid gene; _stop: premature
# stop codon, older versions), sometimes several (e.g. _len_stop). They are not
# part of the sequence ID in the isolate FASTA, so a locus is matched on its raw
# ID first and then with suffixes stripped one at a time (PROGRESS.md R12).
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


def resolve_cell(cell, ids):
    """Resolve each ';'-separated locus of a gene_presence_absence cell."""
    return ";".join(resolve_locus(l, ids) for l in str(cell).split(";"))


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gene_presence_absence",
        help = "CSV file containing binary presence/absence of genes per clone"
    )
    parser.add_argument(
        "--metadata_merged",
        help = "TSV file mapping clone names with sequence data sample names"
    )
    parser.add_argument(
        "--kallisto_dir",
        help = "Directory containing kallisto results",
        default = "./"
    )
    parser.add_argument("--outf", help = "File for results")
    args = parser.parse_args()
    merge_counts(**vars(args))


def merge_counts(gene_presence_absence, metadata_merged, kallisto_dir, outf):
    csv_data = pd.read_csv(gene_presence_absence, low_memory=False)
    metadata = pd.read_csv(metadata_merged, sep="\t")
    colnames = csv_data.columns.values.tolist()
    metadata = metadata[metadata['dna_sample_id'].isin(colnames)]
    quant_dfs = []
    for index, row in metadata.iterrows():
        sample_name = row['sample_name']
        dna_sample_id = row['dna_sample_id']
        quant_file = os.path.join(
            kallisto_dir, 'kallisto_'+sample_name, 'abundance.tsv')
        quant_dat = pd.read_csv(quant_file, sep = "\t")
        quant_dat = quant_dat[["target_id", "est_counts"]]
        quant_dat = quant_dat.rename(columns={'target_id': sample_name})
        cg = csv_data[["Gene", dna_sample_id]]
        cg = cg.rename(columns={dna_sample_id: sample_name})
        ids = set(quant_dat[sample_name])
        cg[sample_name] = cg[sample_name].map(
            lambda c: resolve_cell(c, ids), na_action='ignore')
        all_genes = (cg[sample_name].dropna()).tolist()
        split_genes = [g for g in all_genes if ";" in str(g)]
        split_dict = {} 
        for split_set in split_genes:
            ind_genes = split_set.split(";")
            expr_vals = quant_dat[quant_dat[sample_name].isin(ind_genes)]
            # gene split into several loci in this isolate: total counts over
            # its loci (NaN, not 0, if none of the loci were quantified)
            split_dict[split_set] = expr_vals["est_counts"].sum(min_count=1)
        quant_split = pd.DataFrame([split_dict]).transpose()
        quant_split.rename(columns={0:'est_counts'}, inplace=True)
        quant_split[sample_name] = quant_split.index
        quant_split = quant_split[[sample_name, "est_counts"]]
        quant_split.reset_index(drop = True, inplace = True)
        quant_combined = pd.concat([quant_dat, quant_split], ignore_index = True)
        quant_merged = pd.merge(cg, quant_combined, on=sample_name)
        quant_merged = quant_merged[["Gene", "est_counts"]]
        quant_merged = quant_merged.rename(columns={'est_counts': sample_name})
        quant_dfs.append(quant_merged)
    quant_dfs = [df.set_index('Gene') for df in quant_dfs]
    quant_merged = pd.concat(quant_dfs, axis=1)
    quant_merged.to_csv(outf, index=True, sep='\t')


if __name__ == "__main__":
    parse()

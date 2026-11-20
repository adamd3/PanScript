#!/usr/bin/env Rscript

library(optparse)
library(DESeq2)
if (!require("tidyverse")) {
    lapply(c("readr", "tibble", "purrr", "tidyr", "stringr", "ggplot2"),
        library,
        character.only = TRUE
    )
} else {
    library(tidyverse)
}


option_list <- list(
    make_option(c("-c", "--counts"),
        type = "character", default = NULL,
        help = "table of read counts per gene", metavar = "character"
    ),
    make_option(c("-l", "--lengths"),
        type = "character", default = NULL,
        help = "table of effective lengths per gene", metavar = "character"
    ),
    make_option(c("-m", "--metadata"),
        type = "character", default = NULL,
        help = "sample metadata (sample_name, dna_sample_id); used to take the
                median gene length over isolates, not samples",
        metavar = "character"
    ),
    make_option(c("-o", "--outdir"),
        type = "character", default = NULL,
        help = "output directory for results", metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

counts_f <- opt$counts
lengths_f <- opt$lengths
metadata_f <- opt$metadata
outdir <- opt$outdir
if (is.null(counts_f) || is.null(lengths_f) || is.null(metadata_f) || is.null(outdir)) {
    stop("--counts, --lengths, --metadata and --outdir are all required")
}


## Read data
counts_tab <- suppressMessages(read_tsv(counts_f))
lengths_tab <- suppressMessages(read_tsv(lengths_f))

lengths_tab <- lengths_tab[match(counts_tab$Gene, lengths_tab$Gene), ]

gene_ids <- counts_tab$Gene
## check.names = FALSE keeps sample names as given (e.g. "PA-01" or "1234"
## would otherwise become "PA.01" and "X1234")
counts_tab <- data.frame(counts_tab[, 2:ncol(counts_tab)], check.names = FALSE)
lengths_tab <- data.frame(lengths_tab[, 2:ncol(lengths_tab)], check.names = FALSE)
rownames(counts_tab) <- rownames(lengths_tab) <- gene_ids

## isolate (dna_sample_id) of each sample
metadata <- suppressMessages(read_tsv(metadata_f, col_types = cols(.default = "c")))
isolate <- metadata$dna_sample_id[match(colnames(lengths_tab), metadata$sample_name)]
if (any(is.na(isolate))) {
    stop(paste(
        "samples missing from the metadata:",
        paste(colnames(lengths_tab)[is.na(isolate)], collapse = ", ")
    ))
}

## scale counts to reads per median gene length. The median is over isolates,
## so isolates with several RNA-seq samples don't count more than others:
## first each isolate's length (median over its samples' effective lengths,
## which differ slightly with each sample's fragment length distribution),
## then the median of those over isolates. NA (absent) is ignored throughout.
len_mat <- as.matrix(lengths_tab)
isolate_lens <- matrix(
    vapply(split(seq_len(ncol(len_mat)), isolate),
        function(j) rowMedians(len_mat[, j, drop = FALSE], na.rm = TRUE),
        FUN.VALUE = numeric(nrow(len_mat))
    ),
    nrow = nrow(len_mat)
)
median_lens <- rowMedians(isolate_lens, na.rm = TRUE)
names(median_lens) <- gene_ids

counts_tab_scaled <- (counts_tab / lengths_tab)
counts_tab_scaled <- sweep(counts_tab_scaled, 1, median_lens, "*")


colData <- data.frame(sample_name = colnames(counts_tab_scaled))

# get scaling factors from core gene set
counts_tab_core <- na.omit(counts_tab_scaled)

dds <- suppressMessages(DESeqDataSetFromMatrix(
    countData = round(counts_tab_core),
    colData = colData, design = ~1
))

dds <- estimateSizeFactors(dds)
size_factors <- sizeFactors(dds)



## get size factor-scaled counts (all genes; accessory genes stay NA)
res_df <- sweep(counts_tab_scaled, 2, size_factors, `/`)

## get nTPM values (see docs/normalisation.md):
##   x[g, s] = (counts / isolate's own effective length) / size factor
## i.e. reads per unit length, normalised for library size with the core-gene
## size factors (not per-sample totals, which depend on accessory gene content).
## A single global constant K puts values on a TPM-like scale: per-sample
## totals average 1e6 across samples. K is identical for every gene and sample,
## so it does not affect any comparison.
rate_df <- sweep(counts_tab / lengths_tab, 2, size_factors, `/`)
K <- 1e6 / mean(colSums(rate_df, na.rm = TRUE))
ntpm_df <- rate_df * K


## NB: all outputs are on a linear scale (no log transform); downstream steps
## (e.g. umap.R) apply their own transform where needed.

## convert rownames to column
res_df <- tibble::rownames_to_column(as.data.frame(res_df), "feature_id")
ntpm_df <- tibble::rownames_to_column(as.data.frame(ntpm_df), "feature_id")
counts_tab_scaled <- tibble::rownames_to_column(as.data.frame(
    counts_tab_scaled
), "feature_id")

write.table(
    counts_tab_scaled, file.path(outdir, "raw_counts.tsv"),
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)

write.table(
    res_df, file.path(outdir, "norm_counts.tsv"),
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)

write.table(
    ntpm_df, file.path(outdir, "ntpm.tsv"),
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)

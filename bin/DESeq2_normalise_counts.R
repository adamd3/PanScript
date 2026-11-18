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
    make_option(c("-o", "--outdir"),
        type = "character", default = NULL,
        help = "output directory for results", metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

counts_f <- opt$counts
lengths_f <- opt$lengths
outdir <- opt$outdir


## Read data
counts_tab <- suppressMessages(read_tsv(counts_f))
lengths_tab <- suppressMessages(read_tsv(lengths_f))

lengths_tab <- lengths_tab[match(counts_tab$Gene, lengths_tab$Gene), ]

gene_ids <- counts_tab$Gene
counts_tab <- data.frame(counts_tab[, 2:ncol(counts_tab)])
lengths_tab <- data.frame(lengths_tab[, 2:ncol(lengths_tab)])
rownames(counts_tab) <- rownames(lengths_tab) <- gene_ids

## scale counts to reads per median gene length
median_lens <- rowMedians(as.matrix(lengths_tab), na.rm = TRUE)
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

## get RPKM values
lib_sf <- (colSums(counts_tab_scaled, na.rm = TRUE) / 1e6)
gene_sf <- sapply(median_lens, function(l) {
    l / 1e3 * lib_sf
})
gene_sf <- as.data.frame(t(gene_sf))
rpkm_df <- counts_tab_scaled / gene_sf


## NB: all outputs are on a linear scale (no log transform); downstream steps
## (e.g. umap.R) apply their own transform where needed.

## convert rownames to column
res_df <- tibble::rownames_to_column(as.data.frame(res_df), "feature_id")
rpkm_df <- tibble::rownames_to_column(as.data.frame(rpkm_df), "feature_id")
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
    rpkm_df, file.path(outdir, "rpkm_counts.tsv"),
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)

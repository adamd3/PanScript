process DESEQ_NORMALISE_COUNTS {
    tag "$merged_counts"
    label 'process_medium'
    publishDir "${params.outdir}/gene_counts", mode: 'symlink'

    input:
    tuple path(merged_counts), path(merged_lens)

    output:
    path 'norm_counts.tsv', emit: norm_counts
    path 'ntpm.tsv', emit: ntpm
    path 'raw_counts.tsv', emit: scaled_counts

    script:
    """
    DESeq2_normalise_counts.R \
        -c $merged_counts \
        -l $merged_lens \
        -o ./
    """
}

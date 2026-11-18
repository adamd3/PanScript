process DESEQ_NORMALISE_COUNTS {
    tag "$merged_counts"
    label 'process_medium'
    container 'adamd3/strainseq:latest'
    publishDir "${params.outdir}/gene_counts", mode: 'symlink'

    input:
    tuple path(merged_counts), path(merged_lens)

    output:
    path 'norm_counts.tsv', emit: norm_counts
    path 'rpkm_counts.tsv', emit: rpkm_counts
    path 'raw_counts.tsv', emit: scaled_counts

    script:
    """
    DESeq2_normalise_counts.R \
        -c $merged_counts \
        -l $merged_lens \
        -t TRUE \
        -o ./
    """
}

process UMAP_SAMPLES {
    tag "$norm_counts"
    label 'process_medium'
    publishDir "${params.outdir}/umap_samples", mode: 'copy'

    input:
    path norm_counts
    path meta_merged
    val group
    val seed

    output:
    path '*.{tsv,png}', emit: umap_out

    script:
    """
    umap_samples.py \\
        --norm_counts $norm_counts \\
        --metadata $meta_merged \\
        --group '$group' \\
        --seed $seed \\
        --outdir ./
    """
}

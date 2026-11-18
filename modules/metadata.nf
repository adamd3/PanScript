process CHECK_META_FILE {
    tag "$metadata"
    container 'adamd3/strainseq:latest'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path metadata

    output:
    path 'metadata_final_*.tsv', emit: sample_metadata

    when:
    metadata.exists()

    script:
    def output_file = "metadata_final_${metadata.baseName}.tsv"
    """
    check_metadata.py $metadata ${output_file}
    """
}


/*
 * Report loci listed by Panaroo for an isolate that are missing from its gene
 * FASTA (they cannot be quantified, so their gene groups become NA). Fails if
 * more than params.max_missing_loci of the loci are missing.
 */
process CHECK_FASTA_LOCI {
    tag "$meta.sample_id"
    container 'adamd3/strainseq:latest'

    input:
    path gpa
    tuple val(meta), path(fasta)

    output:
    path "${meta.sample_id}.missing_loci.tsv", emit: missing
    path "${meta.sample_id}.loci_summary.tsv", emit: summary
    stdout emit: message

    script:
    """
    check_fasta_loci.py \\
        --gene_presence_absence $gpa \\
        --dna_sample_id '${meta.dna_sample_id}' \\
        --sample_name '${meta.sample_id}' \\
        --fasta $fasta \\
        --max_missing_frac ${params.max_missing_loci} \\
        --outf ${meta.sample_id}.missing_loci.tsv \\
        --summary ${meta.sample_id}.loci_summary.tsv
    """
}

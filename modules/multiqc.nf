/*
 * One HTML report across all samples: fastp read QC, kallisto pseudo-alignment
 * rates, and the strandedness and missing-loci checks.
 */
process MULTIQC {
    tag "multiqc"
    label 'process_single'
    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    input:
    path inputs, stageAs: 'inputs/?/*'
    path config

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_report_data', emit: data

    script:
    """
    multiqc --force --config $config --filename multiqc_report.html inputs
    """
}

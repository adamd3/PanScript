/*
 * Adaptor trimming and read QC with fastp, using fastp's default settings.
 * Extra options can be passed with --fastp_args (e.g. to approximate the
 * previous Trim Galore defaults: see docs/usage.md).
 */
process FASTP {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/fastp", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith('.fastp.fq.gz')) params.save_trimmed.toString() == 'true' ? "trimmed/$filename" : null
                      else if (filename =~ /\.fastp\.(json|html|log)$/) filename
                      else null   // e.g. the isolate FASTA passed through to quantification
                }

    input:
    tuple val(meta), path(reads), path(fasta)

    output:
    tuple val(meta), path("${meta.sample_id}*.fastp.fq.gz"), path(fasta), emit: trimmed_reads
    path "${meta.sample_id}.fastp.{json,html,log}",                          emit: reports

    script:
    def name    = meta.sample_id
    def args    = params.fastp_args ?: ''
    def threads = Math.min(task.cpus as int, 16)   // fastp uses at most 16 threads
    if (meta.paired_end) {
        """
        fastp \\
            -i ${reads[0]} -I ${reads[1]} \\
            -o ${name}_1.fastp.fq.gz -O ${name}_2.fastp.fq.gz \\
            --thread $threads \\
            --json ${name}.fastp.json --html ${name}.fastp.html \\
            $args \\
            2> ${name}.fastp.log
        """
    } else {
        """
        fastp \\
            -i $reads \\
            -o ${name}.fastp.fq.gz \\
            --thread $threads \\
            --json ${name}.fastp.json --html ${name}.fastp.html \\
            $args \\
            2> ${name}.fastp.log
        """
    }
}

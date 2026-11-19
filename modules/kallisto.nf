process KALLISTO_QUANT {
    tag "$meta.sample_id"
    label 'process_high'
    publishDir "${params.outdir}/kallisto_quant", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads), path(fasta), val(strandedness)

    output:
    path "kallisto_${meta.sample_id}", emit: kallisto_out_dirs

    script:
    def name       = task.ext.prefix ?: "${meta.sample_id}"
    def strand_arg = ['', '--fr-stranded', '--rf-stranded'][strandedness as int]
    // fragment length mean/sd are estimated from paired-end data, but must be
    // given for single-end reads
    def single_arg = meta.paired_end ? '' : "--single -l ${params.fragment_len} -s ${params.fragment_sd}"
    """
    kallisto index -i ${name}.kidx $fasta
    kallisto quant -t $task.cpus -i ${name}.kidx \\
        $strand_arg $single_arg \\
        -o kallisto_${name} $reads
    """
}


process MERGE_COUNTS_AND_LENS {
    tag "merge_counts_and_lens"
    label 'process_single'
    publishDir "${params.outdir}/kallisto_quant", mode: params.publish_dir_mode

    input:
    path gpa_file
    path kallisto_dirs
    path meta_merged

    output:
    tuple path('kallisto_merged_counts_*.tsv'), path('kallisto_merged_lens_*.tsv'), emit: kallisto_merged_out

    script:
    def counts_file = "kallisto_merged_counts_${gpa_file.baseName}_${meta_merged.baseName}.tsv"
    def lens_file = "kallisto_merged_lens_${gpa_file.baseName}_${meta_merged.baseName}.tsv"
    """
    merge_kallisto_counts.py \
        --gene_presence_absence=$gpa_file \
        --metadata_merged=$meta_merged \
        --outf=${counts_file}

    merge_kallisto_lens.py \
        --gene_presence_absence=$gpa_file \
        --metadata_merged=$meta_merged \
        --outf=${lens_file}
    """
}


/*
 * Infer library strandedness from a subsample of reads (first N reads,
 * pseudoaligned with --fr-stranded and with --rf-stranded) and warn if it
 * disagrees with params.strandedness. Does not change any quantification.
 */
process CHECK_STRANDEDNESS {
    tag "$meta.sample_id"
    label 'process_low'

    input:
    tuple val(meta), path(reads), path(fasta)

    output:
    path "${meta.sample_id}.strandedness.tsv", emit: result
    stdout emit: message

    script:
    def name   = meta.sample_id
    def single = meta.paired_end ? "" : "--single -l ${params.fragment_len} -s ${params.fragment_sd}"
    """
    strandedness_check.py subsample --n ${params.strandedness_check_reads} --out sub $reads
    kallisto index -i ${name}.kidx $fasta > kallisto_index.log 2>&1
    for mode in fr rf; do
        kallisto quant -t $task.cpus -i ${name}.kidx --\${mode}-stranded $single \\
            -o q_\${mode} sub_*.fq.gz > kallisto_\${mode}.log 2>&1
    done
    strandedness_check.py infer \\
        --fr q_fr/run_info.json --rf q_rf/run_info.json \\
        --sample '${name}' --specified ${params.strandedness} \\
        --out ${name}.strandedness.tsv
    """
}

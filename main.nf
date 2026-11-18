#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
================================================================================
    StrainQuant: strain-specific analysis of bacterial RNA-Seq data
================================================================================
    Github : [github.com/adamd3/StrainQuant]
*/

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}


/*
================================================================================
    Validate inputs and create channels for files
================================================================================
*/

// Validate required parameters
if (!params.meta_file) {
    exit 1, 'ERROR: --meta_file parameter is required. Please provide a path to the metadata file.'
}
if (!params.gpa_file) {
    exit 1, 'ERROR: --gpa_file parameter is required. Please provide a path to the gene presence/absence file.'
}
if (!params.group) {
    exit 1, 'ERROR: --group parameter is required. Please specify the metadata column for grouping strains.'
}

// Validate parameter values
if (params.strandedness !in [0, 1, 2, 'auto']) {
    exit 1, 'ERROR: --strandedness must be 0 (unstranded), 1 (forward-stranded), 2 (reverse-stranded) or auto (inferred per sample).'
}
if (params.strandedness == 'auto' && params.strandedness_check_reads == 0) {
    exit 1, 'ERROR: --strandedness auto needs the strandedness check; set --strandedness_check_reads > 0.'
}
if (params.fragment_len <= 0) {
    exit 1, 'ERROR: --fragment_len must be greater than 0.'
}
if (params.fragment_sd <= 0) {
    exit 1, 'ERROR: --fragment_sd must be greater than 0.'
}
if (!(params.strandedness_check_reads instanceof Integer) || params.strandedness_check_reads < 0) {
    exit 1, 'ERROR: --strandedness_check_reads must be a whole number >= 0 (0 skips the check).'
}
if (params.max_missing_loci < 0 || params.max_missing_loci > 1) {
    exit 1, 'ERROR: --max_missing_loci must be a fraction between 0 and 1.'
}

// Create file channels with validation
ch_samples = file(params.meta_file, checkIfExists: true)
ch_gpa_file = file(params.gpa_file, checkIfExists: true)




/*
================================================================================
    Modules
================================================================================
*/
include {CHECK_META_FILE; CHECK_FASTA_LOCI} from './modules/metadata'
include {TRIMGALORE} from './modules/trim_reads'
include {KALLISTO_QUANT; MERGE_COUNTS_AND_LENS; CHECK_STRANDEDNESS} from './modules/kallisto'
include {DESEQ_NORMALISE_COUNTS} from './modules/normalisation'
include {UMAP_SAMPLES} from './modules/plots'



/*
================================================================================
    Functions
================================================================================
*/
// Function to get list of [ meta, [ fastq_1, fastq2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create sample metadata
    def meta = [:]
    meta.sample_id     = row.sample_name
    meta.dna_sample_id = row.dna_sample_id
    meta.paired_end    = row.paired.toBoolean()

    // add path(s) of the fastq file(s) to the metadata
    def seq_meta = []
    if (!file(row.fastq1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq1}"
    }
    if (meta.paired_end) {
        if (!file(row.fastq2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq2}"
        }
        seq_meta = [ meta, [ file(row.fastq1), file(row.fastq2) ], file(row.fasta) ]
    } else {
        seq_meta = [ meta, file(row.fastq1), file(row.fasta) ]
    }
    return seq_meta
}


/*
================================================================================
    Main workflow
================================================================================
*/
workflow {

    /*
     *  Create channels for input files
     */

    CHECK_META_FILE (
        ch_samples
    )
    ch_metadata = CHECK_META_FILE.out.sample_metadata

    ch_metadata
        .splitCsv(header: true, sep:'\t')
        .map { create_fastq_channel(it) }
        .set { ch_raw_reads_trimgalore }

    /*
     *  Check that each isolate's Panaroo loci are present in its FASTA
     */
    CHECK_FASTA_LOCI (
        ch_gpa_file,
        ch_raw_reads_trimgalore.map { meta, reads, fasta -> [ meta, fasta ] }
    )
    CHECK_FASTA_LOCI.out.missing
        .collectFile(name: 'missing_loci.tsv', keepHeader: true, skip: 1, sort: true,
                     storeDir: "${params.outdir}/pipeline_info")
    CHECK_FASTA_LOCI.out.summary
        .collectFile(name: 'fasta_loci_summary.tsv', keepHeader: true, skip: 1, sort: true,
                     storeDir: "${params.outdir}/pipeline_info")
    CHECK_FASTA_LOCI.out.message
        .filter { it.trim() }
        .subscribe { log.warn it.trim() }

    /*
     *  Check library strandedness on a subsample of reads (warns on mismatch)
     */
    if (params.strandedness_check_reads > 0) {
        CHECK_STRANDEDNESS (
            ch_raw_reads_trimgalore
        )
        CHECK_STRANDEDNESS.out.result
            .collectFile(name: 'strandedness_check.tsv', keepHeader: true, skip: 1, sort: true,
                         storeDir: "${params.outdir}/pipeline_info")
        CHECK_STRANDEDNESS.out.message
            .filter { it.trim() }
            .subscribe { log.warn it.trim() }
    }

    /*
     *  Trim + pseudo-align reads
     */
    if (params.skip_trimming) {
        ch_reads_for_quant = ch_raw_reads_trimgalore
    } else {
        TRIMGALORE (
            ch_raw_reads_trimgalore
        )
        ch_reads_for_quant = TRIMGALORE.out.trimmed_reads
        ch_trimgalore_results_mqc = TRIMGALORE.out.trimgalore_results_mqc
        ch_trimgalore_fastqc_reports_mqc = TRIMGALORE.out.trimgalore_fastqc_reports_mqc
    }

    // Strandedness for each sample's quantification: the --strandedness value,
    // or (--strandedness auto) the value inferred by CHECK_STRANDEDNESS
    if (params.strandedness == 'auto') {
        ch_inferred_strandedness = CHECK_STRANDEDNESS.out.result
            .map { tsv ->
                def row = tsv.splitCsv(header: true, sep: '\t')[0]
                [ row.sample, row.used_strandedness as int ]
            }
        ch_quant_input = ch_reads_for_quant
            .map { meta, reads, fasta -> [ meta.sample_id, meta, reads, fasta ] }
            .join(ch_inferred_strandedness)
            .map { id, meta, reads, fasta, strandedness -> [ meta, reads, fasta, strandedness ] }
    } else {
        ch_quant_input = ch_reads_for_quant
            .map { meta, reads, fasta -> [ meta, reads, fasta, params.strandedness as int ] }
    }

    KALLISTO_QUANT (
        ch_gpa_file,
        ch_quant_input
    )
    ch_kallisto_out_dirs = KALLISTO_QUANT.out.kallisto_out_dirs

    /*
     *  Merge counts
     */
    MERGE_COUNTS_AND_LENS (
        ch_gpa_file,
        ch_kallisto_out_dirs.collect(),
        ch_metadata
    )
    ch_kallisto_merged_out = MERGE_COUNTS_AND_LENS.out.kallisto_merged_out


    /*
     *  Get size-factor-scaled counts
     */
    DESEQ_NORMALISE_COUNTS (
        ch_kallisto_merged_out
    )
    ch_norm_counts = DESEQ_NORMALISE_COUNTS.out.norm_counts
    ch_scaled_counts = DESEQ_NORMALISE_COUNTS.out.scaled_counts

    /*
     *  UMAP of samples
     */
    UMAP_SAMPLES (
        ch_norm_counts,
        ch_metadata,
        params.group
    )
    ch_umap_out = UMAP_SAMPLES.out.umap_out

}

/*
================================================================================
    Completion summary
================================================================================
*/

c_green = "\033[0;32m";
c_reset = "\033[0m"

workflow.onComplete {
    log.info"""
    Execution status: ${ workflow.success ? 'OK' : 'failed' }
    ${c_green}Results are reported here: $params.outdir${c_reset}
    ￼""".stripIndent()
}

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
      nextflow run StrainQuant --meta_file [file] --gpa_file [gene_presence_absence.csv] --group [str] -profile conda

    Mandatory arguments:
      --meta_file [file]              Path to file containing sample metadata.
      --gpa_file [file]               Path to file containing gene presence/absence per strain (from Panaroo output).
      --group [str]                   Metadata column to use for grouping strains. Default = majority_ST.
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated).
                                      Available: conda, docker

    Other options:
      --strandedness [str]            Strandedness of the reads. 0 = unstranded, 1 = forward-stranded, 2 = reverse-stranded, auto = inferred per sample from a subsample of reads. Default = 2.
      --fragment_len [str]            Estimated average fragment length for kallisto transcript quantification (only required for single-end reads). Default = 150.
      --fragment_sd [str]             Estimated standard deviation of fragment length for kallisto transcript quantification (only required for single-end reads). Default = 20.
      --strandedness_check_reads [int] Reads per sample used to check that --strandedness matches the library (warns on mismatch; results in pipeline_info/strandedness_check.tsv). 0 = skip. Default = 200000.
      --max_missing_loci [float]      Stop if more than this fraction of an isolate's Panaroo loci are missing from its FASTA. Missing loci are always reported in pipeline_info/missing_loci.tsv. Default = 0.5.
      --skip_trimming [bool]          Do not trim adaptors from FastQ files.
      --outdir [file]                 The output directory where the results will be saved (default: './results').
      --cachedir [str]                Defines the path where Singularity images are cached [default: $params.cachedir]
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    """.stripIndent()
}

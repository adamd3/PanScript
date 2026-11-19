#!/usr/bin/env nextflow

/*
================================================================================
    PanScript: isolate-specific bacterial RNA-Seq quantification across a pangenome
================================================================================
    Github : [github.com/adamd3/PanScript]

    Written to parse under both Nextflow 24.10 and the strict syntax of
    Nextflow >= 25 (no top-level statements; errors via error()).
*/

/*
================================================================================
    Modules
================================================================================
*/
include {CHECK_META_FILE; CHECK_FASTA_LOCI} from './modules/metadata'
include {FASTP} from './modules/trim_reads'
include {KALLISTO_QUANT; MERGE_COUNTS_AND_LENS; CHECK_STRANDEDNESS} from './modules/kallisto'
include {DESEQ_NORMALISE_COUNTS} from './modules/normalisation'
include {UMAP_SAMPLES} from './modules/plots'
include {SOFTWARE_VERSIONS} from './modules/versions'



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
        error("ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq1}")
    }
    if (meta.paired_end) {
        if (!file(row.fastq2).exists()) {
            error("ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq2}")
        }
        seq_meta = [ meta, [ file(row.fastq1), file(row.fastq2) ], file(row.fasta) ]
    } else {
        seq_meta = [ meta, file(row.fastq1), file(row.fasta) ]
    }
    return seq_meta
}



// Check parameters before anything runs
def validateParams() {
    if (!params.meta_file) {
        error('ERROR: --meta_file parameter is required. Please provide a path to the metadata file.')
    }
    if (!params.gpa_file) {
        error('ERROR: --gpa_file parameter is required. Please provide a path to the gene presence/absence file.')
    }
    if (!params.group) {
        error('ERROR: --group parameter is required. Please specify the metadata column for grouping strains.')
    }
    if (!(params.seed instanceof Integer)) {
        error('ERROR: --seed must be a whole number.')
    }
    if (params.strandedness !in [0, 1, 2, 'auto']) {
        error('ERROR: --strandedness must be 0 (unstranded), 1 (forward-stranded), 2 (reverse-stranded) or auto (inferred per sample).')
    }
    if (!(params.strandedness_check_reads instanceof Integer) || params.strandedness_check_reads < 0) {
        error('ERROR: --strandedness_check_reads must be a whole number >= 0 (0 skips the check).')
    }
    if (params.strandedness == 'auto' && params.strandedness_check_reads == 0) {
        error('ERROR: --strandedness auto needs the strandedness check; set --strandedness_check_reads > 0.')
    }
    if (params.fragment_len <= 0) {
        error('ERROR: --fragment_len must be greater than 0.')
    }
    if (params.fragment_sd <= 0) {
        error('ERROR: --fragment_sd must be greater than 0.')
    }
    if (params.max_missing_loci < 0 || params.max_missing_loci > 1) {
        error('ERROR: --max_missing_loci must be a fraction between 0 and 1.')
    }
    def publish_modes = ['copy', 'symlink', 'rellink', 'link', 'copyNoFollow', 'move']
    if (!(params.publish_dir_mode in publish_modes)) {
        error("ERROR: --publish_dir_mode must be one of: ${publish_modes.join(', ')}.")
    }
}


/*
================================================================================
    Main workflow
================================================================================
*/
workflow MAIN {

    main:
    ch_samples  = file(params.meta_file, checkIfExists: true)
    ch_gpa_file = file(params.gpa_file, checkIfExists: true)

    SOFTWARE_VERSIONS()


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
        .set { ch_raw_reads }

    /*
     *  Check that each isolate's Panaroo loci are present in its FASTA
     */
    CHECK_FASTA_LOCI (
        ch_gpa_file,
        ch_raw_reads.map { meta, reads, fasta -> [ meta, fasta ] }
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
            ch_raw_reads
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
        ch_reads_for_quant = ch_raw_reads
    } else {
        FASTP (
            ch_raw_reads
        )
        ch_reads_for_quant = FASTP.out.trimmed_reads
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
        params.group,
        params.seed
    )
    ch_umap_out = UMAP_SAMPLES.out.umap_out

}


/*
================================================================================
    Entry point
================================================================================
*/
workflow {
    if (params.help) {
        helpMessage()
    } else {
        validateParams()
        MAIN()
    }
}


def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
      nextflow run adamd3/PanScript --meta_file [file] --gpa_file [gene_presence_absence.csv] -profile docker

    Mandatory arguments:
      --meta_file [file]              Path to file containing sample metadata.
      --gpa_file [file]               Path to file containing gene presence/absence per strain (from Panaroo output).
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated).
                                      Available: docker, singularity, apptainer, conda, test.
                                      Cluster settings: add -c <your_site.config> (see README).

    Other options:
      --group [str]                   Metadata column used to colour samples in the UMAP. Default = majority_ST.
      --strandedness [str]            Strandedness of the reads. 0 = unstranded, 1 = forward-stranded, 2 = reverse-stranded, auto = inferred per sample from a subsample of reads. Default = auto.
      --fragment_len [str]            Estimated average fragment length for kallisto transcript quantification (only required for single-end reads). Default = 150.
      --fragment_sd [str]             Estimated standard deviation of fragment length for kallisto transcript quantification (only required for single-end reads). Default = 20.
      --strandedness_check_reads [int] Reads per sample used to check that --strandedness matches the library (warns on mismatch; results in pipeline_info/strandedness_check.tsv). 0 = skip. Default = 200000.
      --max_missing_loci [float]      Stop if more than this fraction of an isolate's Panaroo loci are missing from its FASTA. Missing loci are always reported in pipeline_info/missing_loci.tsv. Default = 0.5.
      --skip_trimming [bool]          Do not trim reads (skip fastp).
      --fastp_args [str]              Extra options passed to fastp, which otherwise runs with its default settings (e.g. '--cut_tail --cut_tail_mean_quality 20 --length_required 20').
      --save_trimmed [bool]           Save the fastp-trimmed reads in <outdir>/fastp/trimmed.
      --seed [int]                    Random seed for the UMAP, so plots are reproducible. Default = 42.
      --outdir [file]                 The output directory where the results will be saved (default: './results').
      --publish_dir_mode [str]        How results are saved in --outdir: copy, symlink, rellink, link, copyNoFollow or move. Default = copy.
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    """.stripIndent()
}

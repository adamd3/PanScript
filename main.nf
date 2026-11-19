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
include {MULTIQC} from './modules/multiqc'



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



// Parameter values given on the command line are strings in Nextflow >= 25
// (typed in 24.x), so read booleans and numbers through these helpers.
def isTrue(value) {
    return value != null && value.toString().trim().toLowerCase() == 'true'
}

// Convert a parameter value to the first matching JSON-schema type; null if none.
def coerceParam(value, type) {
    def types = type instanceof List ? type : [type]
    def s = value.toString().trim()
    def result = null
    def matched = types.find { t ->
        if (t == 'integer' && (value instanceof Integer || value instanceof Long || s ==~ /[-+]?\d+/)) {
            result = s as Long
        } else if (t == 'number' && (value instanceof Number || s ==~ /[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?/)) {
            result = s as BigDecimal
        } else if (t == 'boolean' && s.toLowerCase() in ['true', 'false']) {
            result = s.toLowerCase() == 'true'
        } else if (t == 'string') {
            result = s
        }
        return result != null
    }
    return matched ? result : null
}

// Validate params against nextflow_schema.json: required, type, allowed values,
// ranges, file existence and unknown (e.g. misspelt) parameters. Done here
// rather than with the nf-schema plugin: the plugin version supporting our
// minimum Nextflow (24.10) rejects command-line numbers on Nextflow >= 25.
def validateParams() {
    def schema = new groovy.json.JsonSlurper().parseText(file("${projectDir}/nextflow_schema.json").text)
    def errors = []
    def known = []
    schema['$defs'].each { key, group ->
        def required = group['required'] ?: []
        group['properties'].each { name, p ->
            known << name
            def raw = params.get(name)
            if (raw == null || raw.toString().trim() == '') {
                if (name in required) {
                    errors << "--${name} is required"
                }
                return
            }
            def v = coerceParam(raw, p['type'])
            if (v == null) {
                errors << "--${name} (${raw}): must be ${p['type'] instanceof List ? p['type'].join(' or ') : p['type']}"
                return
            }
            if (p['enum'] && !p['enum'].any { it.toString() == v.toString() }) {
                errors << "--${name} (${raw}): must be one of ${p['enum'].join(', ')}"
            }
            if (p['minimum'] != null && v < p['minimum']) {
                errors << "--${name} (${raw}): must be at least ${p['minimum']}"
            }
            if (p['maximum'] != null && v > p['maximum']) {
                errors << "--${name} (${raw}): must be at most ${p['maximum']}"
            }
            if (p['exclusiveMinimum'] != null && v <= p['exclusiveMinimum']) {
                errors << "--${name} (${raw}): must be greater than ${p['exclusiveMinimum']}"
            }
            if (p['minLength'] != null && v.toString().length() < p['minLength']) {
                errors << "--${name} must not be empty"
            }
            if (p['exists'] && !file(v.toString()).exists()) {
                errors << "--${name}: file not found: ${v}"
            }
        }
    }
    // Nextflow also stores --a-b style names as a-b and aB
    def variants = known.collectMany { n -> [n, n.replace('_', '-'), n.replaceAll(/_([a-z])/) { all, c -> c.toUpperCase() }] }
    def unknown = params.keySet().findAll { k -> !(k in variants) }
    if (unknown) {
        log.warn "Unrecognised parameter(s), ignored: ${unknown.collect { "--${it}" }.join(', ')} (see --help)"
    }
    // checks across parameters
    if (params.strandedness.toString() == 'auto' && coerceParam(params.strandedness_check_reads, 'integer') == 0) {
        errors << '--strandedness auto needs the strandedness check; set --strandedness_check_reads > 0'
    }
    if (errors) {
        error("Invalid parameters:\n  * " + errors.join("\n  * "))
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
    ch_loci_summary = CHECK_FASTA_LOCI.out.summary
        .collectFile(name: 'fasta_loci_summary.tsv', keepHeader: true, skip: 1, sort: true,
                     storeDir: "${params.outdir}/pipeline_info")
    ch_multiqc_files = Channel.empty().mix(ch_loci_summary)
    CHECK_FASTA_LOCI.out.message
        .filter { it.trim() }
        .subscribe { log.warn it.trim() }

    /*
     *  Check library strandedness on a subsample of reads (warns on mismatch)
     */
    if ((params.strandedness_check_reads as long) > 0) {
        CHECK_STRANDEDNESS (
            ch_raw_reads
        )
        ch_strandedness_tsv = CHECK_STRANDEDNESS.out.result
            .collectFile(name: 'strandedness_check.tsv', keepHeader: true, skip: 1, sort: true,
                         storeDir: "${params.outdir}/pipeline_info")
        ch_multiqc_files = ch_multiqc_files.mix(ch_strandedness_tsv)
        CHECK_STRANDEDNESS.out.message
            .filter { it.trim() }
            .subscribe { log.warn it.trim() }
    }

    /*
     *  Trim + pseudo-align reads
     */
    if (isTrue(params.skip_trimming)) {
        ch_reads_for_quant = ch_raw_reads
    } else {
        FASTP (
            ch_raw_reads
        )
        ch_reads_for_quant = FASTP.out.trimmed_reads
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.reports.flatten().filter { it.name.endsWith('.json') })
    }

    // Strandedness for each sample's quantification: the --strandedness value,
    // or (--strandedness auto) the value inferred by CHECK_STRANDEDNESS
    if (params.strandedness.toString() == 'auto') {
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
    ch_multiqc_files = ch_multiqc_files.mix(KALLISTO_QUANT.out.log)

    /*
     *  QC report across all samples
     */
    if (!isTrue(params.skip_multiqc)) {
        MULTIQC (
            ch_multiqc_files.collect(),
            file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
        )
    }

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
    if (isTrue(params.help)) {
        helpMessage()
    } else {
        validateParams()
        MAIN()
    }
}


// --help, generated from nextflow_schema.json (the single source of parameter
// documentation).
def helpMessage() {
    def schema = new groovy.json.JsonSlurper().parseText(file("${projectDir}/nextflow_schema.json").text)
    def lines = [
        '',
        'Usage:',
        '  nextflow run adamd3/PanScript -profile docker --meta_file metadata.tsv --gpa_file gene_presence_absence.csv --outdir results'
    ]
    schema['$defs'].each { key, group ->
        def required = group['required'] ?: []
        def shown = group['properties'].findAll { name, p -> !p['hidden'] }
        if (shown) {
            lines << ''
            lines << "${group['title']}:"
            shown.each { name, p ->
                def extra = []
                if (name in required) { extra << 'required' }
                if (p['enum']) { extra << "one of: ${p['enum'].join(', ')}" }
                if (p.containsKey('default') && p['default'] != '') { extra << "default: ${p['default']}" }
                def suffix = extra ? " [${extra.join('; ')}]" : ''
                lines << "  --${name.padRight(26)} ${p['description']}${suffix}"
            }
        }
    }
    lines << ''
    lines << 'Software profiles: docker, singularity, apptainer, conda (plus test for the test dataset).'
    lines << 'Cluster settings: add -c <your_site.config> (see docs/usage.md, "Running on a cluster").'
    log.info lines.join('\n')
}

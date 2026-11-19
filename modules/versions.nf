/*
 * Record the versions of the pipeline, Nextflow and every tool used, so a run
 * can be reported (e.g. in a methods section) and reproduced.
 * Written to <outdir>/pipeline_info/software_versions.yml.
 */
process SOFTWARE_VERSIONS {
    tag "software_versions"
    label 'process_low'
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    output:
    path 'software_versions.yml', emit: versions

    script:
    def container = workflow.containerEngine ? (task.container ?: 'none') : 'none (no container engine)'
    """
    cat <<-END_VERSIONS > software_versions.yml
    pipeline:
      name: ${workflow.manifest.name}
      version: ${workflow.manifest.version}
      revision: ${workflow.revision ?: 'local'}
      commit: ${workflow.commitId ?: 'local'}
    nextflow: ${workflow.nextflow.version}
    container_engine: ${workflow.containerEngine ?: 'none'}
    container: ${container}
    tools:
      kallisto: \$(kallisto version | sed 's/.*version //')
      fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
      python: \$(python --version 2>&1 | sed 's/Python //')
      pandas: \$(python -c 'import pandas; print(pandas.__version__)')
      numpy: \$(python -c 'import numpy; print(numpy.__version__)')
      umap-learn: \$(python -c 'from importlib.metadata import version; print(version("umap-learn"))')
      matplotlib: \$(python -c 'import matplotlib; print(matplotlib.__version__)')
      R: \$(Rscript -e 'cat(as.character(getRversion()))')
      DESeq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
    END_VERSIONS
    """
}

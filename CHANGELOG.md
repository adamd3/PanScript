# Changelog

All notable changes to PanScript. Versions follow [semantic versioning](https://semver.org/).

## v2.0.0 (unreleased)

The first release as **PanScript** (previously StrainQuant). Versions before 2026 count as 1.x.

### Changes to outputs

- **New normalised expression measure:** `ntpm.tsv` (normalised TPM) replaces `rpkm_counts.tsv`.
  It's computed with core-genome DESeq2 size factors, so differences in accessory gene content
  between isolates don't distort it (see `docs/normalisation.md`).
- **Linear output matrices:** `norm_counts.tsv` was previously log2(x+1)-transformed; all matrices
  are now on a linear scale. The UMAP log-transforms its input itself.
- **Split genes are summed:** when a gene is split into several loci in an isolate, counts and
  effective lengths are now summed (previously averaged), so the median gene length reflects the
  intact gene. The per-isolate expression rate is unchanged.
- **Panaroo status suffixes are handled:** locus IDs flagged by Panaroo (`_len`, `_pseudo`,
  `_stop`) are now matched to their sequences. Previously these genes were `NA`; on the test data
  this adds 41–47 genes per isolate.
- **All genes are output:** the `--perc` gene filter and the `gene_set_ST_*.tsv` files were
  removed; genes absent from an isolate are `NA`.
- **TMM normalisation removed** (`--norm_method`); DESeq2 size factors only.
- **fastp replaces Trim Galore and FastQC** for trimming and read QC (fastp default settings;
  extra options with `--fastp_args`).
- **UMAP rewritten in Python** (umap-learn) with a fixed `--seed`, so plots are reproducible.
  Coordinates are saved in `umap_coordinates.tsv`.

### New features

- **Strandedness check and `--strandedness auto` (the default):** each library's strandedness is
  inferred from a subsample of reads and used for quantification; a mismatch with an explicitly
  set `--strandedness` gives a warning. Results in `pipeline_info/strandedness_check.tsv`.
- **Missing-loci check:** loci that Panaroo lists but that are missing from an isolate's FASTA
  are reported (`pipeline_info/missing_loci.tsv`, `fasta_loci_summary.tsv`); the run stops if more
  than `--max_missing_loci` of an isolate's loci are missing.
- **MultiQC report** (`multiqc/multiqc_report.html`) across all samples; skip with
  `--skip_multiqc`.
- **Software versions** recorded in `pipeline_info/software_versions.yml`.
- **Parameter validation and `--help`** generated from `nextflow_schema.json`; misspelt
  parameters are reported.
- **`--publish_dir_mode`** (default `copy`): results no longer break if `work/` is deleted.
- **`scripts/prepare_inputs_from_panaroo.py`:** builds the per-isolate gene FASTAs (including
  Panaroo's refound genes) and a metadata template from a Panaroo output directory.
- **Apptainer profile.**

### Breaking changes

- **Nextflow 24.10.0 or later** is required.
- **Single-end and paired-end samples can't be mixed** in one run (the run stops with an error).
- **Site-specific profiles removed** (`cambridge`, `slurm`): put cluster settings in your own
  config and pass it with `-c` (see `docs/usage.md`).
- **`--max_cpus`, `--max_memory` and `--max_time` removed:** use `process.resourceLimits` in a
  config file.
- **`--cachedir` removed:** set `NXF_SINGULARITY_CACHEDIR` / `NXF_APPTAINER_CACHEDIR`.
- **`-profile test` no longer selects Docker:** use `-profile test,docker` (or another software
  profile).

### Reproducibility

- One container image, `adamd3/panscript:2.0.0`, built from `environment.yml` with every tool
  version pinned and always version-tagged (never `latest`).
- CI runs unit tests and the test profile on the minimum supported and latest Nextflow, and
  checks the outputs against a stored baseline.

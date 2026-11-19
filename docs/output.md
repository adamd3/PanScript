# PanScript: output

Everything is written to `--outdir` (default `./results`). By default files are copies
(`--publish_dir_mode copy`), so Nextflow's `work/` directory can be deleted after the run.

- [gene_counts/](#gene_counts): the expression matrices
- [multiqc/](#multiqc): one QC report for the whole run
- [kallisto_quant/](#kallisto_quant): per-sample quantification and merged matrices
- [fastp/](#fastp): read trimming and QC
- [umap_samples/](#umap_samples): sample UMAP
- [pipeline_info/](#pipeline_info): input checks, metadata and software versions

## gene_counts

Three gene × sample matrices, one row per Panaroo gene group and one column per sample
(`sample_name`). All values are on a **linear scale** (not log-transformed). A gene that is absent
from an isolate's genome is `NA`, not 0.

| File | What it is | Compare |
|---|---|---|
| `raw_counts.tsv` | Read counts adjusted **only** for ortholog length differences between isolates (scaled to the gene group's median length). Not normalised for library size. | Nothing directly; it's the input for your own normalisation or DE analysis |
| `norm_counts.tsv` | `raw_counts` divided by DESeq2 size factors estimated on the core genome. | The same gene across samples |
| `ntpm.tsv` | Normalised transcripts per million (nTPM): counts per unit of each isolate's own gene length, divided by the core-genome size factors, scaled so per-sample totals average 10<sup>6</sup>. | Between genes and across samples |

nTPM differs from standard TPM in one way: library size comes from the core-genome size factors,
not from each sample's own total, so differences in accessory gene content between isolates don't
distort it. In practice nTPM is very close to TPM: within a sample it's proportional to kallisto's
TPM, and per-sample totals are approximately rather than exactly 10<sup>6</sup>.

See [normalisation.md](normalisation.md) for the formulas and the reasoning.

## multiqc

`multiqc_report.html`: one report across all samples, made with [MultiQC](https://multiqc.info/):

- **General statistics**: one row per sample with reads kept by fastp, duplication, GC content
  and the kallisto pseudo-alignment rate.
- **fastp**: read quality, adapter content and filtering, before and after trimming.
- **kallisto**: reads processed and pseudo-aligned per sample. A low rate can mean the wrong
  strandedness, contamination or the wrong FASTA.
- **Strandedness check** and **Panaroo loci missing from FASTA**: the input-check tables from
  `pipeline_info/` (below).

`multiqc_report_data/` has the same numbers as text files. Skip the report with `--skip_multiqc`.

## kallisto_quant

- `kallisto_<sample>/`: kallisto's output for each sample (`abundance.tsv` with `est_counts`,
  `eff_length` and `tpm` per locus of that isolate, and `run_info.json`).
- `<sample>.kallisto.log`: kallisto's log for each sample (reads processed and pseudo-aligned).
- `kallisto_merged_counts_*.tsv` and `kallisto_merged_lens_*.tsv`: estimated counts and effective
  lengths mapped from loci to Panaroo gene groups (gene group × sample). When a gene is split into
  several loci in an isolate, its counts and lengths are summed. These are the inputs to the
  normalisation; most users want `gene_counts/` instead.

## fastp

`<sample>.fastp.html`, `.json` and `.log`: [fastp](https://github.com/OpenGene/fastp)'s report
for each sample. With `--save_trimmed`, the trimmed reads are saved in `fastp/trimmed/`. Not
produced with `--skip_trimming`.

## umap_samples

- `umap_<group>.png`: UMAP of the samples, coloured by the `--group` metadata column.
- `umap_coordinates.tsv`: the coordinates behind the plot.

The UMAP is computed on log2(`norm_counts` + 1) of the core genes (present in every sample), with
`n_neighbors` = 15 (or 3 with fewer than 15 samples, and always fewer than the number of samples),
and a fixed `--seed`, so it's reproducible.

## pipeline_info

| File | What it is |
|---|---|
| `metadata_final_*.tsv` | The metadata as used by the run (with a `paired` column added). |
| `missing_loci.tsv` | One row per locus that Panaroo lists for an isolate but that isn't in its FASTA (empty if none). Those loci can't be quantified. Flags: `refound_stop` / `refound_len` (Panaroo "refound" genes), `split_partial` (only some loci of a split gene are missing, so it's undercounted), `split_all_missing`. |
| `fasta_loci_summary.tsv` | Per sample: loci listed, loci missing and gene groups affected. |
| `strandedness_check.tsv` | Per sample: reads tested, the share pseudo-aligning in each orientation, the inferred strandedness, the strandedness specified, whether they agree, and the strandedness used for quantification. |
| `software_versions.yml` | Versions of the pipeline, Nextflow, the container image and every tool used. Useful for a methods section. |

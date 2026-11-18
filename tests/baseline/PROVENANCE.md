# Baseline provenance

These outputs record pipeline behaviour **before** the 2026 publication changes (PROGRESS.md T0).

| Item | Value |
|------|-------|
| Date | 2026-11-18 |
| Pipeline code | identical to `main` @ `87be5d8` (run from `dev` @ `f186581`, which differs only in documentation) |
| Command | `nextflow run . -profile test,docker --outdir <results> -w <work> -with-trace` |
| Nextflow | 23.10.1 (Java 11, macOS arm64; container under amd64 emulation) |
| Container | `adamd3/strainseq@sha256:0e549c29a0a139749d069ed0411601650fd33036140ad3aec939b3a3b318dd12` (amd64, built 2025-07-29) |
| Tools in container | kallisto 0.51.1 · Trim Galore 0.6.10 · cutadapt 5.1 · FastQC 0.12.1 · Python 3.11.13 · pandas 2.3.1 · R 4.4.3 · DESeq2 1.46.0 · edgeR 4.4.0 · umap 0.2.10.0 · matrixStats 1.5.0 |
| Inputs | `configs/test.config`, which downloads `test_data/` from GitHub `main` (3 *P. aeruginosa* isolates, single-end, about 100k reads each) |
| Parameters | test-profile defaults: `strandedness=2`, `perc=99`, `norm_method=DESeq`, `fragment_len=150`, `fragment_sd=20`, `group=majority_ST` |

## Contents
- `gene_counts/`: `raw_counts.tsv`, `norm_counts.tsv` (log2(x+1)), `rpkm_counts.tsv` (log2(x+1)),
  and `gene_set_ST_*` (the unused `--perc` subset)
- `kallisto_quant/`: merged `est_counts` and `eff_length` matrices, plus each sample's
  `abundance.tsv` and `run_info.json`
- `pipeline_info/`: processed metadata and the Nextflow trace
- `umap_samples/`: UMAP object and plot

## Known properties of this baseline (see PROGRESS.md R9, R10)
- **Only 2.6–6.1% of reads pseudo-align.** The test library is forward-stranded, but the default
  `strandedness=2` (reverse) is used. With `strandedness=1` the rate is about 50%. So the counts
  are very sparse, and DESeq2 size factors come out at 0.928 / 1.000 / 1.000.
- 41–47 genes per isolate are NA (treated as absent) even though Panaroo lists them as present.
  These are `*_refound_*_stop` loci, which aren't in the input FASTA files.

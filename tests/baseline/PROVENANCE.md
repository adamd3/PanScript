# Baseline provenance

These outputs record pipeline behaviour **before** the 2026 publication changes (PROGRESS.md T0).

| Item | Value |
|------|-------|
| Date | 2026-11-18 |
| Pipeline code | identical to `main` @ `87be5d8`, except that `configs/test.config` sets `strandedness = 1` (see History) |
| Command | `nextflow run . -profile test,docker --outdir <results> -w <work> -with-trace` |
| Nextflow | 23.10.1 (Java 11, macOS arm64; container under amd64 emulation) |
| Container | `adamd3/strainseq@sha256:0e549c29a0a139749d069ed0411601650fd33036140ad3aec939b3a3b318dd12` (amd64, built 2025-07-29) |
| Tools in container | kallisto 0.51.1 · Trim Galore 0.6.10 · cutadapt 5.1 · FastQC 0.12.1 · Python 3.11.13 · pandas 2.3.1 · R 4.4.3 · DESeq2 1.46.0 · edgeR 4.4.0 · umap 0.2.10.0 · matrixStats 1.5.0 |
| Inputs | `configs/test.config`, which downloads `test_data/` from GitHub `main` (3 *P. aeruginosa* isolates, single-end, about 100k reads each) |
| Parameters | `strandedness=1` (forward), `perc=99`, `norm_method=DESeq`, `fragment_len=150`, `fragment_sd=20`, `group=majority_ST` |

## Contents
- `gene_counts/`: `raw_counts.tsv`, `norm_counts.tsv` (log2(x+1)), `rpkm_counts.tsv` (log2(x+1)),
  and `gene_set_ST_*` (the unused `--perc` subset)
- `kallisto_quant/`: merged `est_counts` and `eff_length` matrices, plus each sample's
  `abundance.tsv` and `run_info.json`
- `pipeline_info/`: processed metadata and the Nextflow trace
- `umap_samples/`: UMAP object and plot

## Known properties of this baseline
- Pseudo-alignment: 51.0% (ZG205864), 50.1% (ZG301975), 46.9% (ZG302367).
- **The size factors (1.002 / 1.000 / 1.000) can't be trusted as a check on normalisation.** At
  about 100k reads per sample, the median rounded count per core gene is 1. 56% of core genes
  contain a zero and DESeq2 excludes them, and about 5% of the remaining ratios are exactly 1
  (genes with equal counts in all samples), which is where the median lands. Normalisation maths
  is therefore checked with synthetic data in `tests/unit/`.
- 41–47 genes per isolate are NA (treated as absent) even though Panaroo lists them as present.
  These are `*_refound_*_stop` loci, which aren't in the input FASTA files (PROGRESS.md R9/P7).

## History
- **v0** (commit `cf2084e`): same code, but the test profile used the default `strandedness=2`.
  Only 2.6–6.1% of reads pseudo-aligned, because the test library is forward-stranded
  (PROGRESS.md R10). Replaced on 2026-11-18 by the current baseline once the test profile was
  set to `strandedness = 1`. Retrieve it with `git show cf2084e:tests/baseline/<file>`.

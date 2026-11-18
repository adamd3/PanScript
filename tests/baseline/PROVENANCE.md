# Baseline provenance

The reference outputs that later changes are checked against (see `tests/README.md`).

**Current baseline: v3, after Phase 2 (fastp, Python UMAP, pinned container) with
`--strandedness auto` as the default.**

| Item | Value |
|------|-------|
| Date | 2026-11-18 |
| Pipeline code | `feat/phase2-container` (merged into `dev` on 2026-11-18) |
| Command | `nextflow run . -profile test,docker --outdir <results> -w <work>` (no overrides) |
| Nextflow | 23.10.1 (Java 11, macOS arm64; container under amd64 emulation) |
| Container | `adamd3/panscript:2.0.0-dev` = `sha256:b97d5bec60caf57a5a68ca6a16151dbc76d3d511a8847abd4c6d588d8f2e3fb3` (identical image first published as `adamd3/strainseq:2.0.0-dev` before the rename to PanScript) (amd64; built from this repo's `environment.yml` + `Dockerfile`) |
| Tools in container | kallisto 0.51.1 · fastp 1.4.0 · Python 3.11.13 · pandas 2.3.1 · numpy 2.3.2 · umap-learn 0.6.0 · matplotlib 3.10.9 · R 4.4.3 · DESeq2 1.46.0 · optparse 1.7.5 · matrixStats 1.5.0 · tidyverse 2.0.0 |
| Inputs | `configs/test.config`, which downloads `test_data/` from GitHub `main` (3 *P. aeruginosa* isolates, single-end, about 100k reads each) |
| Parameters | defaults: `strandedness=auto` (inferred forward for all 3 samples), `fragment_len=150`, `fragment_sd=20`, `group=majority_ST`, `seed=42`, fastp defaults |

## Contents
- `gene_counts/`: `raw_counts.tsv`, `norm_counts.tsv`, `ntpm.tsv` (all linear; NA = absent)
- `kallisto_quant/`: merged `est_counts` and `eff_length` matrices, plus each sample's
  `abundance.tsv` and `run_info.json`
- `fastp/`: fastp JSON reports
- `pipeline_info/`: processed metadata, `missing_loci.tsv`, `fasta_loci_summary.tsv`,
  `strandedness_check.tsv`
- `umap_samples/`: `umap_coordinates.tsv` (deterministic with the fixed seed, so it's compared)
  and the plot

## Known properties of this baseline
- fastp keeps 99.0–99.4% of reads. Pseudo-alignment: 51.1% (ZG205864), 50.7% (ZG301975),
  47.4% (ZG302367).
- **The size factors can't be trusted as a check on normalisation.** At about 100k reads per
  sample, most core genes have counts of 0–2. Normalisation maths is checked with synthetic data in
  `tests/unit/`.
- 41–47 genes per isolate are NA because their Panaroo loci (`*_refound_*_stop` / `_len`) aren't in
  the input FASTA files. They're listed in `pipeline_info/missing_loci.tsv`.

## History
Retrieve an old version with `git show <commit>:tests/baseline/<file>`.

- **v0** (commit `cf2084e`): code identical to `main` @ `87be5d8`, test profile with
  `strandedness=2`. Only 2.6–6.1% of reads pseudo-aligned, because the library is forward-stranded.
- **v1** (commit `7ff012e`): code identical to `main` @ `87be5d8`, with `strandedness = 1`. This is the
  pre-publication-changes reference (Trim Galore, log2 outputs, averaged split genes, RPKM).
- **v2** (commit `a067dd2`): after the Phase 1 normalisation changes, using the old
  `adamd3/strainseq:latest` image (`sha256:0e549c29…`) with Trim Galore 0.6.10.
- **v3** (this version): after Phase 2. Compared with v2: fastp instead of Trim Galore (merged-count
  Spearman 0.999, nTPM 0.997; same genes and NA pattern), the UMAP is now Python with a fixed seed,
  and all other tool versions are unchanged. `--strandedness auto` gives the same result as v2's
  explicit `1`.

# Baseline provenance

The reference outputs that later changes are checked against (see `tests/README.md`).

**Current baseline: v2, after the Phase 1 normalisation changes (PROGRESS.md P2–P6).**

| Item | Value |
|------|-------|
| Date | 2026-11-18 |
| Pipeline code | `feat/normalisation` @ `24704ab` |
| Command | `nextflow run . -profile test,docker --outdir <results> -w <work>` |
| Nextflow | 23.10.1 (Java 11, macOS arm64; container under amd64 emulation) |
| Container | `adamd3/strainseq@sha256:0e549c29a0a139749d069ed0411601650fd33036140ad3aec939b3a3b318dd12` (amd64, built 2025-07-29) |
| Tools in container | kallisto 0.51.1 · Trim Galore 0.6.10 · cutadapt 5.1 · FastQC 0.12.1 · Python 3.11.13 · pandas 2.3.1 · R 4.4.3 · DESeq2 1.46.0 · umap 0.2.10.0 · matrixStats 1.5.0 |
| Inputs | `configs/test.config`, which downloads `test_data/` from GitHub `main` (3 *P. aeruginosa* isolates, single-end, about 100k reads each) |
| Parameters | `strandedness=1` (forward), `fragment_len=150`, `fragment_sd=20`, `group=majority_ST` |

## Contents
- `gene_counts/`: `raw_counts.tsv`, `norm_counts.tsv`, `ntpm.tsv` (all linear; NA = absent)
- `kallisto_quant/`: merged `est_counts` and `eff_length` matrices (split loci summed), plus each
  sample's `abundance.tsv` and `run_info.json`
- `pipeline_info/`: processed metadata
- `umap_samples/`: UMAP object and plot. **Not deterministic between runs** (PROGRESS.md R11), so
  don't compare these.

## Known properties of this baseline
- Pseudo-alignment: 51.0% (ZG205864), 50.1% (ZG301975), 46.9% (ZG302367).
- **The size factors (≈1.002 / 1.000 / 1.000) can't be trusted as a check on normalisation.** At
  about 100k reads per sample, the median rounded count per core gene is 1. 56% of core genes
  contain a zero and DESeq2 excludes them, and the median ratio lands on exact ties. Normalisation
  maths is checked with synthetic data in `tests/unit/`.
- 41–47 genes per isolate are NA (treated as absent) even though Panaroo lists them as present.
  These are `*_refound_*_stop` loci, which aren't in the input FASTA files (PROGRESS.md R9/P7).

## History
Retrieve an old version with `git show <commit>:tests/baseline/<file>`.

- **v0** (commit `cf2084e`): code identical to `main` @ `87be5d8`, with the test profile's default
  `strandedness=2`. Only 2.6–6.1% of reads pseudo-aligned, because the library is forward-stranded
  (R10).
- **v1** (commit `7ff012e`): code identical to `main` @ `87be5d8`, with `strandedness = 1`. This is the
  pre-publication-changes reference. Outputs included `gene_set_ST_*.tsv`, log2(x+1)
  `norm_counts`/`rpkm_counts`, and averaged split genes. Every Phase 1 step was checked against v1;
  the observed differences are recorded in PROGRESS.md (P2–P6).
- **v2** (this version): replaced v1 on 2026-11-18 once Phase 1 was complete. It matches the
  final Phase 1 run exactly.

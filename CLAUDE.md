# CLAUDE.md

Instructions for AI agents (and humans) working on this repo. Read this file first, then read
[PROGRESS.md](PROGRESS.md) before starting any work.

## What this project is

A Nextflow (DSL2) pipeline that quantifies gene expression in bacterial RNA-Seq data **without a
single reference genome**. Each input "isolate" (a genomically distinct unit: a strain, lineage,
clone or isolate; we don't assume they're strains) has its own gene set. Reads are pseudo-aligned
with kallisto to that isolate's own genes. Counts are then merged across isolates through
orthologous gene groups from a Panaroo pangenome (`gene_presence_absence` file) and normalised.

Name: **PanScript** (decided 2026-11-18; previously StrainQuant, originally OrthoQuant). Docker
image: `adamd3/panscript` (lowercase, as Docker requires). Target: publication in **Genome Biology**.

### Pipeline flow (main.nf)
1. `CHECK_META_FILE`: validates metadata and adds a `paired` column (`bin/check_metadata.py`); then
   `CHECK_FASTA_LOCI` reports Panaroo loci missing from each isolate FASTA (`bin/check_fasta_loci.py`)
   and `CHECK_STRANDEDNESS` warns if `--strandedness` disagrees with the reads (`bin/strandedness_check.py`)
2. `FASTP` (skipped with `--skip_trimming`): adaptor trimming and read QC, fastp defaults plus `--fastp_args`
3. `KALLISTO_QUANT`: builds a per-isolate index from that isolate's `fasta` and quantifies
4. `MERGE_COUNTS_AND_LENS`: maps kallisto targets to Panaroo gene groups, giving merged
   `est_counts` and `eff_length` matrices (`bin/merge_kallisto_{counts,lens}.py`)
5. `DESEQ_NORMALISE_COUNTS`: length-scaling, size factors, nTPM
   (`bin/DESeq2_normalise_counts.R`)
6. `UMAP_SAMPLES`: UMAP of samples coloured by `--group`, fixed `--seed` (`bin/umap_samples.py`)

User-run helper scripts (not part of the workflow) live in `scripts/`:
`prepare_inputs_from_panaroo.py` (per-isolate FASTA + metadata template from a Panaroo run) and
`build_and_push_image.sh` (versioned image release).

## Ground rules

1. **Consult the user before any change that affects the output of any pipeline step.** That
   includes numeric values, file names, columns, formats, default parameters and filtering. You may
   suggest alternatives to existing methods, but only when you think they'll genuinely help. Explain
   why, then wait for approval. Refactors, docs, CI and tooling that don't change outputs are fine to
   do directly.
2. **Don't assume unusual-looking code is a mistake.** Many choices are deliberate (see below).
   If you're unsure why something is done, ask or log it as an open question in PROGRESS.md.
   Don't "fix" it silently.
3. **Keep PROGRESS.md up to date.** Several agents may work in parallel. Claim a task before you
   start it, update it when you finish, and record decisions there (see the protocol at the top of
   PROGRESS.md).
4. Only commit or push when the user asks. Follow the git workflow below.
5. **Never add AI attribution to commits or PRs.** That means no `Co-Authored-By: Claude` trailer,
   no "Generated with Claude Code" line and no session links. Commit messages describe the change
   only.
6. Write in British English in user-facing text (normalise, colour), matching the existing code.

## Git workflow

- `main` is the stable, released pipeline. It's what `nextflow run adamd3/<repo>` pulls by
  default. Don't commit to it directly.
- `dev` is the integration branch. Branch each task off `dev` as `feat/<short-name>`,
  `fix/<short-name>` or `docs/<short-name>`, then merge back into `dev`. **Until CI is fixed
  (PROGRESS.md H3), merge locally rather than through a GitHub PR.** Every PR currently triggers
  a rebuild and push of the image, which changes the tools used by every run. (No longer true since H3: CI never pushes images, so PRs are safe.) Give each task its own
  branch, so parallel agents don't collide. Note the branch in PROGRESS.md.
- To release, merge `dev` into `main` and tag it `vX.Y.Z` (semantic versioning). Tag the Docker
  image the same way.
- (This dev/master split is also nf-core's convention, so it costs nothing if we go nf-core later.)

## Deliberate design choices (don't change these without discussion)

- **Ortholog length scaling.** An orthologous group can differ in length between isolates. Counts
  are divided by each isolate's kallisto effective length and multiplied by the median length of
  that gene across isolates. This gives "reads per median gene length".
- **Core-gene size factors.** Accessory genes are NA in the merged matrix. Size factors (DESeq2
  median-of-ratios) are estimated only on genes with no NAs, which are the core genes
  present in every sample. They're then applied to all non-NA values in each sample. This assumes
  the core genome (usually thousands of genes) is big enough to estimate library size accurately.
- **NA means absent, not zero.** A gene missing from an isolate's genome is NA, not a zero count.
- **Split genes.** When a Panaroo cell holds several `;`-separated loci (a gene split into 2+
  genes in some isolates but not others), their counts and effective lengths are both **summed**
  over the quantified loci (since 2026-11-18; previously averaged). The value is NA if none of
  the loci were quantified.
- **Per-sample kallisto index.** Each sample is indexed against its own isolate's gene FASTA.
- **DESeq2 only** (decided 2026-11-18; TMM removed). The `--perc` gene filter was also
  removed: the pipeline outputs every gene, with NA for genes absent from an isolate.

### Meaning of the output matrices (agreed with user 2026-11-18)
The full rationale and formulas are in [docs/normalisation.md](docs/normalisation.md). Keep that
file in sync with any change to normalisation code.

- `raw_counts.tsv`: counts adjusted **only** for ortholog length varying between isolates (scaled
  to median length). Not normalised for library size, and not comparable between genes.
- `norm_counts.tsv`: `raw_counts` divided by DESeq2 size factors from core genes. Adjusted for
  ortholog length variation and library size. Comparable across samples within a gene, but not
  between different genes.
- `ntpm.tsv` (nTPM; replaced `rpkm_counts.tsv`): (count / the isolate's own effective length) /
  size factor, times one global constant to put it on a TPM-like scale. Comparable both between
  genes and across samples.
- All matrices are **linear** (no log transform).
- **UMAP log-transforms its input itself** (log2(x+1) in `bin/umap_samples.py`). Don't feed it already
  logged values, and don't remove that line.

## Repo practicalities

- Notes for one particular machine (Java and Docker quirks) go in `CLAUDE.local.md`, which is
  gitignored, not here.

- Test run: `nextflow run . -profile test,docker --outdir ./results`. Test inputs are fetched from
  GitHub `main` via raw URLs (`configs/test.config`, `test_data/meta_merged.txt`), so CI tests the
  data on `main`, not the data on your branch. The GitHub repo is `adamd3/PanScript`
  (renamed from `adamd3/StrainQuant`; the old names `StrainQuant`, `OrthoQuant` and `StrainSeq`
  redirect to it).
- In a run, **`bin/` scripts come from your local checkout** (Nextflow puts `bin/` on PATH inside
  the container), but **tools (kallisto, fastp, R and its packages, umap-learn) come from the Docker
  image**. So script changes are tested straight away, but a new tool (for example fastp) needs a
  rebuilt image or a per-process container.
- Most `publishDir`s use `mode: 'symlink'`, so `results/` files point into `work/`. Use `cp -L`
  when saving outputs, and don't delete `work/` before copying.
- CI (`.github/workflows/ci.yml`, on pushes to `main`/`dev` and on PRs) runs the unit tests inside
  the pinned image, then the test profile on Nextflow 23.10.1 and latest-stable, and checks the
  outputs against `tests/baseline` with `tests/compare_to_baseline.py`. **CI never pushes
  images.** `container.yml` only checks that the image still builds when `environment.yml` or the
  `Dockerfile` changes. A PR that changes outputs on purpose must re-capture the baseline in the
  same PR.
- **Changing the software environment:** edit `environment.yml`; run
  `scripts/build_and_push_image.sh <version>`, which builds, smoke-tests and runs the unit tests in
  the image, then asks before pushing; **ask the user** before confirming the push; then set that tag in `nextflow.config`
  (`process.container`) and re-capture the baseline.
- **Nextflow versions:** 23.10.0 and later (`manifest.nextflowVersion`). Code must parse under
  Nextflow's strict syntax (default from 25.x): no statements outside workflows or functions; use
  `error()`, not `exit`; no helper functions in config files; no `$VAR` environment variables in
  config strings.
- Containers: every process uses the single `process.container` set in `nextflow.config` (no
  per-module containers), built from `environment.yml` with **exact pinned versions** and the
  micromamba-based `Dockerfile`. Run the unit tests inside the image (`docker run --rm -v
  $PWD:/repo -w /repo <image> python -m unittest discover -s tests/unit`), because some tests
  need umap-learn.
- Despite its `.tsv` extension, `test_data/gene_presence_absence.tsv` is a comma-separated Panaroo
  CSV. The scripts read it with `pd.read_csv` (comma-delimited).
- **Panaroo locus IDs:** in `gene_presence_absence.csv`, Panaroo appends status suffixes (`_len`,
  `_pseudo`, or `_stop` in older versions) to locus IDs. These are not part of the sequence IDs in
  `gene_data.csv` or the isolate FASTA, so strip them when matching (see PROGRESS.md R12).
- `KALLISTO_QUANT` takes reads as given (trimmed by fastp, or raw with `--skip_trimming`) and
  strandedness per sample in its input tuple (fixed value, or inferred with `--strandedness auto`).

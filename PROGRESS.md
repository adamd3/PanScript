# PROGRESS.md

Shared progress tracker for everyone working on this repo, people and agents alike.
Read [CLAUDE.md](CLAUDE.md) first.

## Protocol

- **Claim a task before starting.** Change its status to `IN PROGRESS (<agent/branch>, YYYY-MM-DD)`.
  Don't start a task someone else has claimed. If two tasks touch the same files (for example
  P2/P3/P4 all edit the normalisation R scripts), do them on one branch, or one after another.
- **When you finish**, set the status to `DONE (YYYY-MM-DD)`, add a one-line summary, and name the
  branch, PR or commit.
- **Decisions** made with the user go in the Decision log. Write down the reasoning, not just the
  result.
- **Questions for the user** go in Open questions. Don't act on an output-affecting change until
  the user has answered.
- Status values: `TODO` · `IN PROGRESS` · `BLOCKED (reason)` · `NEEDS USER (question)` · `DONE`.
- Use absolute dates. Branch each task off `dev` (see Git workflow in CLAUDE.md).

---

## ▶ Start here (status as of 2026-11-19, end of day)

**To resume:** read CLAUDE.md, then this section. Work is on `dev` (clean, pushed, CI green).
No task is claimed or in progress, and there are no unmerged branches.

**v2.0.0 is released.** Tag `v2.0.0` on `main`; GitHub release
https://github.com/adamd3/PanScript/releases/tag/v2.0.0; image `adamd3/panscript:2.0.0`
(`sha256:7a1fbe6c…`); Zenodo concept DOI **10.5281/zenodo.24375403** (v2.0.0:
10.5281/zenodo.24375404). Verified: `nextflow run adamd3/PanScript -r v2.0.0 -profile test,docker`
from GitHub matches baseline v4. Phases 0–4 are done.

**Since the release, on `dev` only:** the DOI badge and citation in README, CITATION.cff and
CITATIONS.md (commit 078b454). `main` doesn't show the DOI yet (see Q15).

**User's machine:** fixed on 2026-11-19. Plain `nextflow` is now Homebrew 26.10.6 (conda's
Nextflow was removed and `JAVA_HOME` in `~/.bash_profile` was corrected).

**Next steps, in order:**
1. **Ask the user Q15** (merge the DOI docs commit into `main` now, or wait for the next release?).
2. **Manuscript (Phase 5):** start with M1 (outline for a Genome Biology Software/Method
   article: check current article types and requirements, and do a literature search on related
   tools). Then M2 (benchmark design, including PanXpress and microberna comparisons; needs the
   user's input), M3 (showcase dataset), M4 (redraw the pipeline diagram;
   `docs/images/panscript_pipeline_lowres.png` predates fastp/MultiQC/input checks).
3. **CITATION.cff co-authors** (the user will supply them later).
4. **Smaller open items:** R7 (document the recommended downstream DE workflow), H3 remainder (a
   test option using local `test_data/` instead of GitHub URLs).

For later releases: bump `manifest.version`, build/push a new image with
`scripts/build_and_push_image.sh`, update `process.container`, CHANGELOG, merge `dev` → `main`,
tag, and `gh release create` (Zenodo mints a new version DOI automatically).

---

## Plan (overview)

Status: Phase 0 ✓ · Phase 1 ✓ · Phase 1b ✓ · Phase 2 ✓ · Phase 3 ✓ (except the H3 local-test option) ·
Phase 4 ✓ (v2.0.0 released 2026-11-19) · Phase 5 (manuscript) not started.

Phases run in order. Within a phase, tasks on separate branches can run in parallel.

**Phase 0: Baseline (T0).** Before any code changes, run the test profile on unchanged `main` and
save its outputs, as a fixed record of current behaviour. Then every later change can be checked
against it.

**Phase 1: Normalisation changes.** These are all on one branch, `feat/normalisation`, because
they edit the same scripts. Do them in this order, checking against the baseline after each one:
1. P3: remove `--perc` and `SUBSET_GENES`. Expected: outputs identical (the subset was never used).
2. P4: remove TMM. Expected: DESeq2 outputs identical.
3. P6: linear outputs, with `umap.R` applying log2 itself. Expected: new = 2^old − 1, and the UMAP
   is identical.
4. P5: sum split genes. Expected: only split-gene rows change in `raw_counts`/`norm_counts`, and
   their c/L rate is unchanged.
5. P2: replace RPKM with nTPM. New output; check by hand on synthetic data.
6. Update `docs/normalisation.md` (remove the *[planned]* markers) and the README Output section.

**Phase 1b: Input checks.** P7 (missing-loci report) and P8 (strandedness check). These add new
reports and warnings but don't change existing values, so they can go on their own branches in
parallel with Phase 1.

**Phase 2: Container rebuild bundle, P1 (fastp) + P9 (Python UMAP) + H2 (pinned versions).** One rebuild, one baseline re-capture. Originally: trimming (P1, fastp), on `feat/fastp`. This can start in parallel with Phase 1. It
needs fastp in the container: rebuild the image or use a biocontainer per process. Expected:
counts change slightly (different trimmer), so compare correlations with the baseline rather than
requiring identical values.

**Phase 3: Reproducibility and repo hygiene.** H2 (pin versions and version-tag the image), H3
(fix CI, plus a test option that reads local files), H4/H7/H8 (cleanup), H5 (README), then
H6 (nf-core cost/benefit write-up).

**Phase 4: Name and release (H1).** Rename everything (done: PanScript). Merge `dev` into
`main` and tag `v2.0.0` (the manifest is at 2.0.0dev; the pre-2026 code counts as 1.x). Get a Zenodo
DOI.

**Phase 5: Manuscript (M1–M4).** Outline and benchmarking design can start any time, but the
results need the Phase 1–2 pipeline.

---

## Tasks

### Pipeline changes (these affect outputs, so get user sign-off on the details)

| # | Task | Status | Notes |
|---|------|--------|-------|
| P1 | Replace Trim Galore with fastp for adaptor trimming; drop FastQC (fastp does QC). | DONE (2026-11-18, `feat/phase2-container`, merged into `dev`) | **Result:** `FASTP` replaces Trim Galore and FastQC (fastp 1.4.0 defaults, `--fastp_args` passthrough, reports in `<outdir>/fastp/`, `--save_trimmed` → `fastp/trimmed/`). KALLISTO_QUANT no longer depends on Trim Galore file names (and its unused `gpa` input was removed). Test data vs baseline v2: fastp keeps 99.0–99.4% of reads; pseudoalignment +0.1 to +0.6 points; same genes and NA pattern; Spearman 0.999 (merged counts, norm_counts) and 0.997 (nTPM); median gene change 0–1%; 1.2–1.7% of expressed genes change by > 10% (mostly low counts). **Decided 2026-11-18:** run fastp with **its own defaults**, which are more intuitive for new users. Add a `--fastp_args` passthrough for extra options, and document in the README how to approximate Trim Galore (for example `--cut_tail --cut_tail_mean_quality 20 --length_required 20`). Note: fastp's defaults include read-level quality filtering (reads with > 40% of bases below Q15 are discarded), a 15 bp minimum length and automatic polyG trimming for NovaSeq/NextSeq data. Approved. Also update the `KALLISTO_QUANT` filename symlinks (they rely on Trim Galore's names), `environment.yml`, README and `save_trimmed` publishing. Consider publishing the fastp JSON/HTML reports, and optionally MultiQC (already in env). Agree fastp parameters with the user. |
| T0 | Capture baseline outputs and add a regression check. | DONE (2026-11-18, `feat/regression-baseline`, merged into `dev`) | The baseline is in `tests/baseline/` (provenance in `PROVENANCE.md`), and `tests/compare_to_baseline.py` checks new runs against it (it matches the baseline itself exactly). Confirmed R1, R2, R4 and the log2 outputs on real data. Found **R9** and **R10**. Still to do: synthetic-data unit tests (T1). |
| T1 | Unit tests on synthetic data for the merge and normalisation scripts. | DONE (2026-11-18, `feat/normalisation`) | **Result:** 14 tests in `tests/unit/test_normalisation.py` (merging, split genes including unquantified/partial loci, NA handling, size factors, norm_counts, nTPM formula and properties, and an accessory-gene robustness test showing why nTPM rather than TPM). Checked that they fail on the old averaging code. `tests/unit/test_normalisation.py`: 7 tests, all passing against *current* behaviour, using an independent numpy reference. Update the expected values with each Phase 1 step. | The test dataset is too sparse to check normalisation maths (see R10). Build tiny hand-computable inputs (split genes, accessory NA, known size factors) and check each output formula, especially nTPM. Write it before or during Phase 1. |
| P2 | Replace RPKM with nTPM, computed from DESeq2 core-gene size factors. | DONE (2026-11-18, `feat/normalisation`) | **Result:** `ntpm.tsv` replaces `rpkm_counts.tsv`. Test data: per-sample totals 0.94–1.06 × 1e6 (mean exactly 1e6); within each sample nTPM ÷ kallisto TPM is constant (CV ~1e-6) at 0.94–1.07, so it is TPM up to a per-sample size-factor adjustment. `raw_counts`/`norm_counts` unchanged by this step. Name still to be confirmed (Q1b). Approved 2026-11-18 (option b). Definition: x_gs = (c_gs / L_gs) / sf_s, where L is the isolate's own effective length, then multiplied by one global constant K so that the mean over samples of the per-sample total (all non-NA genes) is 1e6. K is identical for every gene and sample, so it doesn't affect comparisons. Proposed name: **nTPM** ("normalised TPM", a precedent from the Human Protein Atlas), file `ntpm.tsv`; check the name with the user. The README and methods must state that per-sample totals are only *approximately* 1e6. Median length cancels out, so P5 doesn't affect this output. |
| P3 | Remove `--perc`, `SUBSET_GENES` and `gene_set_ST_*.tsv`; output every gene. | DONE (2026-11-18, `feat/normalisation`) | **Result:** `gene_set_ST_*.tsv` no longer produced; all other outputs identical to the baseline. Approved 2026-11-18. Also remove the `-p` branch in the R scripts, `bin/subset_genes.py`, and the param from config, help, README and validation. Accessory genes stay NA. |
| P4 | Remove TMM: DESeq2 only. | DONE (2026-11-18, `feat/normalisation`) | **Result:** test outputs identical to the baseline. Approved 2026-11-18. Remove `TMM_NORMALISE_COUNTS`, `bin/TMM_normalise_counts.R` and `--norm_method` (help, config, README, validation). edgeR can then come out of `environment.yml`. |
| P5 | Sum (not average) counts and effective lengths for split genes. | DONE (2026-11-18, `feat/normalisation`) | **Result:** merged counts and lengths = baseline × number of loci for split cells only; 87 `raw_counts` rows changed, all split genes (median length typically ×2); rate c/L unchanged for every gene. DESeq2's rounding of counts before estimating size factors moved one size factor by 0.02%, so `norm_counts` shifts by that factor in one column. Approved 2026-11-18. Edit `bin/merge_kallisto_counts.py` and `bin/merge_kallisto_lens.py`. The per-isolate rate is unchanged; the median length, and so the magnitude of `raw_counts` and `norm_counts`, changes only for genes that are split in some isolates. |
| P6 | Output linear (non-log) matrices. | DONE (2026-11-18, `feat/normalisation`) | **Result:** new = 2^baseline − 1 for norm/rpkm; UMAP input identical (max |diff| 1e-14). UMAP itself is not deterministic between runs (see R11). Approved 2026-11-18. Drop the hard-coded `-t TRUE`, and say clearly in the README (and in file names or headers if useful) that values are untransformed. **`UMAP_SAMPLES` currently relies on log2 `norm_counts`, so `umap.R` must apply log2(x+1) itself** to keep the plot unchanged. |
| P7 | Report loci that Panaroo lists but that are missing from the isolate FASTA. | DONE (2026-11-18, `feat/missing-loci`) | **Result:** `bin/check_fasta_loci.py` + `CHECK_FASTA_LOCI` (per sample, alongside trimming). Writes `pipeline_info/missing_loci.tsv` (flags `refound_stop`, `refound_len`, `split_partial`, `split_all_missing`) and `pipeline_info/fasta_loci_summary.tsv`, and prints a `WARN` line per affected sample. New param `--max_missing_loci` (default 0.5) stops the run above that fraction; tested in a real run. Test data: 41/42/47 loci missing (0.6–0.7%), all refound `_stop` (129) or `_len` (1). Expression outputs identical to baseline v2. 6 unit tests. Approved 2026-11-18 (R9/Q6): use whatever sequences the user provides, and **report** what's missing. Don't change values (missing loci stay NA, as now). Proposed design: a per-sample check before kallisto compares the Panaroo loci for that `dna_sample_id` with the FASTA headers and writes `pipeline_info/missing_loci.tsv` (sample, gene group, locus, flags such as `refound_stop` and `partial_split`, meaning only some loci of a split cell are present), plus a warning with counts per sample. Proposal (needs user OK, since it can stop runs): **fail** if most loci are missing (for example > 50%), which almost certainly means the wrong FASTA or a mismatched `dna_sample_id`. Note that with P5 (summing), a `partial_split` gene undercounts, so it must be flagged. |
| P8 | Automatic strandedness check. | DONE (2026-11-18, `feat/strandedness-check`) | **Result:** `bin/strandedness_check.py` + `CHECK_STRANDEDNESS` (per sample, on the first `--strandedness_check_reads` reads, default 200k, 0 = skip; runs alongside trimming). Thresholds follow nf-core/rnaseq: forward if ≥ 80% of stranded pseudoalignments are forward, reverse if ≤ 20%, unstranded if 40–60%, else undetermined. Writes `pipeline_info/strandedness_check.tsv`, and WARNs on mismatch with the suggested setting. **Warn-only: never changes quantification** (`--strandedness auto` not built → Q10). Test data: forward (89–95%), agrees with `strandedness = 1`; with `--strandedness 2` all three samples warn, as expected. Expression outputs identical to baseline v2. 9 unit tests (they caught a floating-point boundary bug at exactly 20%). Approved 2026-11-18 (Q7b). Quantify a subsample of reads (for example 200k) per sample in both orientations (kallisto `--fr-stranded` / `--rf-stranded`). Infer forward/reverse/unstranded from the ratio, and **warn** when it disagrees with `--strandedness` (and possibly offer `--strandedness auto` to use the inferred value, as nf-core/rnaseq does). Report the result per sample in `pipeline_info/`. Thresholds and the auto mode need user sign-off. |
| P9 | UMAP: fixed random seed, rewritten in Python. | DONE (2026-11-18, `feat/phase2-container`, merged into `dev`) | **Result:** `bin/umap_samples.py` (umap-learn 0.6.0, `compatibility_layout=True`, `n_jobs=1`, `--seed` default 42) replaces `bin/umap.R`. Outputs `umap_coordinates.tsv` + `umap_<group>.png`. Identical coordinates across runs (verified). With 3 samples `n_neighbors` is capped at 2 (umap-learn needs k < n). Numeric labels tidied ('313.0' → '313'). 3 unit tests (run in the container). Approved 2026-11-18 (Q9). Replace `bin/umap.R` with a Python script using `umap-learn` with a fixed `random_state` (new `--seed` param, default 42), plus matplotlib or seaborn for the plot. Keep the same input (log2(`norm_counts` + 1), core genes only, NA rows dropped), the same `n_neighbors` rule (3 if fewer than 15 samples, else 15), and the same outputs (a saved embedding, which can no longer be `.rds`, so use a TSV of coordinates, plus `umap_<group>.png`). Needs `umap-learn` and a plotting library in the container. The plot will differ from the R version (different implementation and seed), which is expected. Check: two runs give identical coordinates. |
| P10 | `--strandedness auto`. | DONE (2026-11-18, `feat/strandedness-auto`) | **Result:** `--strandedness auto` uses `CHECK_STRANDEDNESS`'s inferred value per sample (joined on sample name into KALLISTO_QUANT, which now takes strandedness in its input tuple), and stops with an error if undetermined. New column `used_strandedness` in `strandedness_check.tsv`. Default is still `2` (Q11 open). Test runs: default → identical to baseline v2; `auto` → inferred forward for all 3 samples, outputs identical to baseline v2; `--skip_trimming` path still works. 3 new unit tests (29 total). Approved 2026-11-18 (Q10). Use the strandedness inferred by `CHECK_STRANDEDNESS` for each sample's main kallisto run. Stop with a clear error if it's undetermined. Record the value used in `pipeline_info/strandedness_check.tsv`. KALLISTO_QUANT must then wait for the check (join on sample). Open question Q11: should `auto` become the default? |
| — | Suggested order: T0 → one branch `feat/normalisation` for P3 + P4 + P5 + P6 + P2 (they all touch the same R and Python scripts) → P1 (fastp; independent, can run in parallel on its own branch) → H-tasks. | | |

### Publication / repo hygiene (no change to outputs)

| # | Task | Status | Notes |
|---|------|--------|-------|
| H1 | Choose a new pipeline name. | DONE (2026-11-18, `chore/rename-panscript`, merged into `dev`); GitHub repo renamed to `adamd3/PanScript` by the user, local remote updated | **Decided 2026-11-18: PanScript.** Renamed everywhere describing the current pipeline (README, main.nf, manifest, environment name, test config and test-data URLs, CLAUDE.md, diagram filename); historical records keep old names. Image re-tagged and pushed as `adamd3/panscript:2.0.0-dev` (same digest `b97d5bec…`). URLs point to `github.com/adamd3/PanScript` (verified 200; old names redirect). Docs diagram `panscript_pipeline_lowres.png` still shows the old name and Trim Galore steps → redraw (M4). History: OrthoQuant (original) → StrainQuant. The GitHub repo is currently `adamd3/StrainQuant`, and the old names `OrthoQuant` and `StrainSeq` redirect to it. Renaming the repo back to OrthoQuant would work smoothly with GitHub's redirects. Candidates proposed 2026-11-18: **OrthoQuant** (recommended), **PanQuant**. A quick web search found no tool with either exact name; the nearest are *OrthoQuantum* (eukaryotic phylogenetic profiling, NAR 2022) and *PanXpress* (bacterial pan-transcriptomics). Before deciding, check bioconda, PyPI and GitHub for clashes. It should describe what the pipeline does: per-isolate pseudo-alignment and pangenome-aware merging of expression. It shouldn't imply "strains". Then rename everything: repo, README, `main.nf` header and help, Docker image `adamd3/strainseq`, conda env name, CI, test URLs. |
| H2 | Pin software versions and rebuild the container. | MOSTLY DONE (2026-11-18) | **Published** `adamd3/strainseq:2.0.0-dev`, now `adamd3/panscript:2.0.0-dev` (`sha256:b97d5bec…`) with the user's approval; it's set once in `nextflow.config` and `--pull=always` was removed. Baseline v3 was captured with it. **Remaining:** a final release tag (with the final name, H1), an optional multi-arch build, and CI still builds and pushes `:latest` (H3). **Done:** `environment.yml` fully pinned (kallisto, DESeq2, R and pandas unchanged from the baseline image; adds fastp 1.4.0, umap-learn 0.6.0, matplotlib 3.10.9, procps-ng; drops `defaults`, Trim Galore, FastQC, edgeR, R umap, multiqc and pyfaidx); micromamba 1.5.10 Dockerfile without `COPY . /app`; `.dockerignore`; per-module containers removed (one `process.container` in nextflow.config). Built locally as `strainquant-local:phase2` (amd64); all 32 unit tests pass inside it. **Remaining:** publish the image with a version tag (needs the user's OK and the final name), set that tag in config, drop `--pull=always`, re-capture baseline v3 with the published digest, multi-arch build (optional). **Not before T0; combine the rebuild with P1** (fastp has to go in the image anyway). Dockerfile/env issues: (1) base `continuumio/miniconda3` is unpinned; pin it or switch to a pinned `mambaorg/micromamba`. (2) `environment.yml` uses `>=`; pin exact versions. (3) Drop the Anaconda `defaults` channel (licensing terms for large organisations); use conda-forge + bioconda with strict priority. (4) `COPY . /app` bakes the whole repo (test data, stray files) into the image; it isn't needed, since Nextflow supplies `bin/`. Remove it or add a `.dockerignore`. (5) Update dependencies: add fastp; remove trim-galore, fastqc and (after P4) edgeR; keep multiqc only if used. (6) Replace `:latest` and `--pull=always` (always-pull also makes every task contact Docker Hub: on 2026-11-18 a test run failed on a registry timeout and needed `-resume`) with a version tag `adamd3/<name>:vX.Y.Z`, set once in config (modules also hard-code the image). (7) Build for amd64 **and** arm64 with `docker buildx` *if* all bioconda dependencies have linux-aarch64 builds; check each one, otherwise stay amd64-only. After rebuilding, compare outputs with the T0 baseline: a new kallisto or DESeq2 version could shift values slightly. |
| H3 | Fix CI. | DONE (2026-11-18, `fix/ci`); **first CI run on `dev` green**: unit tests, plus 23.10.1 and latest-stable both matching baseline v3 exactly on Linux, UMAP included | **Result:** new `ci.yml`: unit tests inside the pinned image; test profile on Nextflow 23.10.1 and latest-stable; outputs compared against `tests/baseline`. Triggers: push to main/dev, PRs. No secrets, **no image pushes**. New `container.yml`: build-only check when environment.yml/Dockerfile change. **Nextflow 26 compatibility (found while testing locally):** the strict syntax rejected `check_max()` in config (replaced with inline closures), `$PWD` in config (`System.getenv('PWD')`), top-level statements in main.nf (moved into `validateParams()` + `workflow MAIN`), `exit` (`error()`), and `workflow.onComplete` in the body (removed; Nextflow prints its own summary). `when: metadata.exists()` made Nextflow 26 silently skip **every** task (URL input), so it was removed. `docker.userEmulation` is gone in 26, so `docker.runOptions = '-u $(id -u):$(id -g)'` keeps tasks running as the calling user (the new image's default user is `mambauser`, which would break on Linux). Removed `nextflow.enable.dsl=2` (minimum is now 23.10, where DSL2 is the default). Added a `manifest` (nextflowVersion >= 23.10.0, version 2.0.0dev). Verified locally: 23.10.1 and 26.10.6 both give outputs identical to baseline v3; `--help` and parameter errors work on both. Still open: test data comes from GitHub `main` URLs; optional local test option. No run history shows on GitHub. That's expected: Actions is enabled and the secrets exist, but the last push was in July 2025, so old runs have passed the retention limit. The push job's `if:` checks `github.repository == 'adamd3/StrainQuant'`. That's correct today, but it must be updated if the repo is renamed (H1). **Side effect to fix before using PRs:** every PR run rebuilds the image from unpinned `environment.yml` and **pushes it as `adamd3/strainseq:latest`**. That silently changes the tools that every `-profile docker` run (and the T0 baseline digest) uses. Until H3 is fixed, merge branches into `dev` locally instead of through GitHub PRs. It pushes `:latest` on every PR, and fork PRs have no secrets. It only triggers on pushes to `main` (add `dev`). It uses old actions (`checkout@v2`, `login-action@v1`). Test data comes from `main`, not the branch under test. The minimum Nextflow version is 21.10.3, which is old. Consider adding small unit tests (pytest/testthat) for the merging and normalisation scripts. |
| H4b | Test `scripts/prepare_inputs_from_panaroo.py` on a real Panaroo run. | DONE (2026-11-19, `feat/panaroo-helper-gz`) | **Result:** user's run (2,967 genomes; `gene_data.csv.gz` 10.7 GB, 19.8M genes of which 1.66M refound). For the 3 test isolates the helper's FASTAs are **identical** to `test_data/SRR87372{81,82,83}.fna` (IDs and sequences), and every locus in `gene_presence_absence.csv` has a sequence. Added gzip input support; gpa loci are now kept only for the selected isolates (peak memory 3.1 GB → 43 MB; 137 s). 2 new unit tests. `gene_data.csv` stays the source: it's the only standard Panaroo output with each isolate's own sequence for refound genes (`pan_genome_reference.fa` has one sequence per group; `combined_DNA_CDS.fasta` is, as far as we know, written before refinding and uses internal clustering IDs; not verified against Panaroo's source). Limitation: without `--isolates`, all sequences are held in memory before writing. The big files are excluded locally via `.git/info/exclude`, never commit them. |
| H4 | Clean up unused code. | DONE (2026-11-18, `chore/cleanup`) | **Result:** removed `modules/make_clone_fasta.nf`, `bin/make_single_clone_fasta.py`, `bin/merge_metadata.py`, `bin/length_scale_counts.R` (its median-length scaling already lives in `DESeq2_normalise_counts.R`), `configs/container.config` and `requirements.txt` (user approved). **Replaced** the clone-FASTA helper with `scripts/prepare_inputs_from_panaroo.py`, which uses Panaroo's `gene_data.csv` (includes refound genes), strips Panaroo status suffixes and writes a metadata template. 4 unit tests on a synthetic Panaroo directory built from Panaroo's own writer code (not yet tested on a real Panaroo run, since no `gene_data.csv` is in test_data). See "Unused files" in CLAUDE.md. Confirm with the user before deleting (some may be useful as helper utilities). |
| H5 | Update README. | TODO | Use the agreed definitions of the output matrices (CLAUDE.md, "Meaning of the output matrices"). Update parameters after P3/P4. Add a citation section. Do this after P1–P6 land. |
| H6 | Assess nf-core compliance. | DONE (2026-11-19, `docs/nf-core-assessment`): write-up in `docs/nf-core-assessment.md`. Recommendation: no listing before submission; take the cheap wins (schema + nf-schema, software versions file, `--publish_dir_mode copy`, MultiQC at the next image rebuild, docs layout). Cheap wins approved (Q13) → H9–H14. | The user is reluctant to do heavy lifting just for nf-core. First write up a cost/benefit: what nf-core listing needs (template, nf-schema, modules, MultiQC, nf-test) versus cheap wins worth taking anyway (`nextflow_schema.json`, `-profile test`, versioned containers, `versions.yml`). |
| H15 | Minimum Nextflow 24.10; `resourceLimits`; right-sized resources (Q14). | DONE (2026-11-19, `feat/resources`) | **Result:** `manifest.nextflowVersion = '!>=24.10.0'` (23.10 now stops with a clear message); `max_cpus`/`max_memory`/`max_time` params and their closures replaced by `process.resourceLimits` (test profile: 2 CPUs, 6 GB). Labels: default 1 CPU/4 GB; new `process_single` 1 CPU/8 GB (CHECK_FASTA_LOCI, MERGE, DESeq2, UMAP: all single-threaded); `process_low` 2/4 GB (strandedness check); `process_medium` 4/8 GB (fastp); `process_high` 8/8 GB (kallisto; was 12 CPUs/72 GB). Memory doubles on an OOM retry. Also: `check_fasta_loci.py` and the merge scripts now read only the needed columns of gene_presence_absence. On the user's 2,967-genome table: check 4.8 GB → 1.7 GB (21 s → 4 s), merge 3.9 GB → 1.8 GB (17 s → 4 s), identical output. CI matrix now 24.10.0 + latest-stable. Verified: 24.10.0 and 26.10.6 both match baseline v4; 45/45 unit tests in the image. |
| H9 | `nextflow_schema.json`: parameter validation and generated `--help`. | DONE (2026-11-19, `feat/schema`) | **Result:** `nextflow_schema.json` (JSON schema 2020-12, nf-core layout, so Seqera Platform can build a launch form) documents every parameter. `validateParams()` in `main.nf` validates against it (required, type, enum, ranges, file exists, unrecognised/misspelt params → warning) and `--help` is generated from it. **The nf-schema plugin was tried and dropped:** 2.4.2 (newest supporting Nextflow 24.10) ignores its config on Nextflow 26 and rejects numeric command-line values there; unpinned, Nextflow 26 gets 3.1.0 but every run prints a large "pin your version" warning. Our own ~60-line validator works identically on 24.10.0 and 26.10.6 and needs no plugin download (good for offline clusters). **Found and fixed a pre-existing bug:** on Nextflow ≥ 25 all command-line params are strings (and `params` can't be reassigned), so e.g. `--skip_trimming false` was truthy and skipped trimming. Booleans/numbers are now read via `isTrue()` / `as long` / `.toString()`. Verified: 24.10.0 plain and 26.10.6 with every default given on the command line both match baseline v4; `--skip_trimming` on 26 skips fastp. |
| H10 | Record tool versions in `pipeline_info/software_versions.yml`. | DONE (2026-11-19, `feat/publish-copy`): new `SOFTWARE_VERSIONS` process (`modules/versions.nf`) writes pipeline name/version/commit, Nextflow version, container engine and image, and kallisto, fastp, Python, pandas, numpy, umap-learn, matplotlib, R, DESeq2 versions. | Approved (Q13). New output file. kallisto, fastp, Python, pandas, numpy, umap-learn, R, DESeq2. |
| H11 | `--publish_dir_mode`, default `copy`. | DONE (2026-11-19, `feat/publish-copy`): every `publishDir` uses it; validated in `validateParams()`; README, help and CLAUDE.md updated. Matrices identical to baseline v4; no symlinks in `results/`. | Approved (Q13). Outputs become copies, not symlinks into `work/`. Update CLAUDE.md's `cp -L` note and the README. |
| H12 | MultiQC report (fastp + kallisto). | DONE (2026-11-19, `release/v2.0.0`) | **Result:** `MULTIQC` process (`modules/multiqc.nf`, `assets/multiqc_config.yml`), MultiQC 1.35 in `environment.yml`. Sections: general stats, fastp, kallisto (from new `kallisto_quant/<sample>.kallisto.log`), and our strandedness and missing-loci tables as custom content. Sample names aligned across sections. `--skip_multiqc`. Tested with the MultiQC biocontainer before the image rebuild; expression matrices unchanged. |
| H13 | nf-core-style docs. | DONE (2026-11-19, `release/v2.0.0`) | **Result:** README slimmed to intro, quick start, outputs summary and citation; details moved to `docs/usage.md` and new `docs/output.md`; `CITATIONS.md`, `CHANGELOG.md` (v2.0.0 vs 1.x, incl. breaking changes) and `CITATION.cff` added. **NEEDS USER:** ORCID / co-authors for CITATION.cff. |
| H14 | Config overhaul: no site-specific profiles. | DONE (2026-11-19, `feat/config-overhaul`) | **Result:** removed the `cambridge` and `slurm` profiles (user: relics; site settings go in a `-c` config, documented in README "Running on a cluster"), the unused `name`, `cachedir` and `config_profile_*` params, and the Singularity `--bind $PWD`/cacheDir (use `NXF_SINGULARITY_CACHEDIR`). Added an `apptainer` profile. The `test` profile no longer picks the container engine itself: use `-profile test,docker` (or singularity/apptainer/conda). Outputs identical to baseline v4. Resource limits (`max_*` params vs `resourceLimits`) and right-sizing the process labels wait on Q14. |
| H7 | Clean up the repo. | DONE (2026-11-19: CITATION.cff, Zenodo DOI at the v2.0.0 release) — was MOSTLY DONE (2026-11-18, `chore/cleanup`) | **Result:** deleted the 3 empty `test_data/.!NNNNN!*.fq.gz` files; replaced the untracked `docker_build_and_push.sh` (pushed `:latest`) with `scripts/build_and_push_image.sh <version>` (reads the repo from nextflow.config, builds, smoke-tests, runs unit tests in the image, asks before pushing). **Remaining:** CITATION.cff and Zenodo DOI (at release). Stray zero-byte `test_data/.!NNNNN!*.fq.gz` files (failed-transfer leftovers, untracked). Decide whether to commit `docker_build_and_push.sh`. Add a CITATION.cff and a Zenodo DOI before submission. |
| H8 | Fix minor code issues (R8). | TODO | No output change. |

### Manuscript (Genome Biology)

| # | Task | Status | Notes |
|---|------|--------|-------|
| M1 | Draft an outline for a Genome Biology "Method" or "Software" article. | TODO | Check the current GB article types and requirements: availability, OSI licence (MIT ✓), archived version. |
| M2 | Plan benchmarking and validation. | TODO | **Related tools to compare against or discuss:** PanXpress (PyPI `panxpress`: pan-transcriptomic gapped k-mer index, alignment-free, includes its own annotation-correction step); the Greene lab "2022-microberna" single-pangenome versus best-genome comparison (github.com/greenelab/2022-microberna); single-reference mapping. Run a full literature search during M1. For example, compare against single-reference mapping (such as PAO1) on multi-isolate datasets; use simulated reads with known expression for different isolate sets and accessory fractions; check robustness of core-gene size factors as core size shrinks. The design needs input from the user. |
| M3 | Choose the showcase dataset(s). | TODO | Test data is *P. aeruginosa* (SRX5123741/2/4, the HZI AMR cohort). |
| M4 | Redraw the pipeline diagram (`docs/images/panscript_pipeline_lowres.png`). | TODO | The current image shows the old name and the old steps (Trim Galore/FastQC, RPKM). It should show input checks, fastp, per-isolate kallisto, ortholog merging, core-genome size factors, nTPM, UMAP. Useful as manuscript Figure 1. |

### T0 procedure (baseline run)

How the test is wired (checked 2026-11-18): `-profile test` loads `configs/test.config`. That
sets `meta_file` and `gpa_file` to raw GitHub URLs for `test_data/` on `main`.
`meta_merged.txt` in turn points `fastq1`/`fasta` at raw GitHub URLs. So the test uses the
`test_data/` files, but **downloads them from GitHub `main`**, not from the local folder. That's the
same input CI uses: 3 *P. aeruginosa* isolates, single-end, subsampled. Local `main` matches
`origin/main` and all the URLs return 200.

Steps:
1. Run on unchanged `main` (code untouched; only the notes are committed):
   `nextflow run . -profile test,docker --outdir <scratch>/baseline_results -w <scratch>/work`
   (with `JAVA_HOME` set as described in CLAUDE.local.md on the user's Mac).
2. Record the Nextflow version, git commit and the image digest of `adamd3/strainseq:latest`
   (amd64 only; last pushed 2025-07-29). `--pull=always` means a later image push could change
   the tools, which is why the digest matters.
3. Copy the outputs with `cp -L` (they're symlinks into `work/`) to `tests/baseline/`. The key
   files are `gene_counts/{raw_counts,norm_counts,rpkm_counts}.tsv`, the
   `kallisto_quant/kallisto_merged_*` tables and `pipeline_info/metadata_final_*.tsv`. They're
   small enough to commit.
4. Write a small comparison script (`tests/compare_to_baseline.py`) that checks each Phase 1
   expectation and reports any unexpected differences.

Later (H3): add a test option that reads the **local** `test_data/` files, so changes to test
data on a branch get tested, and runs need no network.

---

## Code review findings (2026-11-18)

- **R1. `--perc` does nothing.** Resolved: remove it (P3).
- **R2. RPKM ignores size factors.** It uses column sums of length-scaled counts over all non-NA
  genes. Resolved: being replaced (P2).
- **R3. TMM `norm_counts` are on a different scale from DESeq2's.** Resolved: TMM is being removed
  (P4).
- **R4. Split genes.** Counts and effective lengths are both averaged across the `;`-separated
  loci, so counts/length = sum/sum and the per-isolate rate is correct. **But** the median length
  across isolates is taken from the *averaged* (per-fragment) length. If a gene is split into two
  in most isolates, the median comes out at about half the full gene length. This doesn't affect
  comparisons across samples within a gene (it's a per-gene constant), but it does make
  `raw_counts` too small for that gene and would bias any between-gene measure (P2). Fix: sum
  instead of averaging (P5). Small additional caveat: kallisto's effective length subtracts about
  one fragment length per locus, so two fragments have a summed effective length slightly shorter
  than the intact gene. This is minor and probably not worth correcting.
- **R5. Paired-end detection is global.** `check_metadata.py` marks every sample as paired if any
  `fastq2` is filled in, so mixed single-end and paired-end inputs would break. Unanswered: should
  mixed input be supported, or rejected with a clear error? **[DONE 2026-11-19, `fix/reject-mixed-reads`]**
  Rejected (Q4): `check_metadata.py` stops with an error naming the single-end and paired-end
  samples. 3 unit tests (`tests/unit/test_check_metadata.py`); test metadata output unchanged.
- **R6. Log2 outputs.** `norm_counts.tsv` and `rpkm_counts.tsv` are log2(x+1)-transformed
  (`-t TRUE` is hard-coded), but the file names and README don't say so (P6).
- **R7. Downstream differential expression.** Length scaling gives non-integer `raw_counts`. If
  users take them into DESeq2 for differential expression, document the recommended route:
  rounded `raw_counts`, or raw `est_counts` plus a per-gene, per-sample length offset (as tximport
  does). Unanswered; lower priority, but reviewers may ask.
- **R9. [Resolved by P7: now reported] Pseudogenised "refound" genes become NA.** Found 2026-11-18 in the baseline. Panaroo
  lists `*_refound_*_stop` loci (refound genes with a premature stop codon) as present, but the
  input isolate FASTA files leave them out. Refound loci *without* `_stop` are included. So
  kallisto never quantifies them, and they end up NA, which is treated as "absent". In the test
  data: 129 such loci, 41–47 genes per isolate, and 39 genes present in all 3 isolates drop out
  of the core set used for size factors. Leaving them out may be deliberate (pseudogenes). But
  NA conflates "absent from genome" with "present but pseudogenised". The FASTA files are built
  outside the pipeline, so this depends on how users make them. → Q6.
- **R10. Test data is forward-stranded, but the default is `--strandedness 2` (reverse).** Found
  2026-11-18. Re-running kallisto on the test reads gives pseudo-alignment of about 51% (forward),
  about 55% (unstranded) and **3–6% (reverse, the default)**. So the test (and CI) has always run
  on about 10% of its reads, making the baseline counts very sparse. The default (reverse, as in
  dUTP/TruSeq stranded kits) is sensible in general, but it's wrong for this dataset. **If any
  real analyses of this cohort used the default, they lost about 90% of reads.** → Q7.
- **R11. UMAP is not reproducible between runs.** Found 2026-11-18. `umap.R` uses `umap.defaults`, whose `random_state` is drawn at random. Two runs on identical input gave layouts differing by up to 0.97 units. For a publication, figures should be reproducible: set a fixed `random_state` (perhaps as a `--seed` parameter). This changes the plot, so it needs the user's OK → Q9.
- **R12. "Missing" loci are really an ID-matching bug.** Found 2026-11-18 while reading Panaroo's source (`generate_output.py`). In `gene_presence_absence.csv`, Panaroo appends status suffixes to locus IDs (`_len` if the length deviates from the group mode, `_pseudo` for invalid genes, `_stop` in older versions). These suffixes are **not part of the sequence ID** (`gene_presence_absence_roary.csv` has the plain IDs). In the test data, **all 130 "missing" loci (R9/P7) are in the isolate FASTA** under the plain ID, for example `1_refound_1076_stop` is the sequence `1_refound_1076`. So 41–47 genes per isolate (39 of them core) are NA only because the merge scripts compare IDs literally. Fix: strip the suffix when matching (raw ID first, then stripped), in `merge_kallisto_*.py` and `check_fasta_loci.py`. **Changes outputs → Q12.** **[DONE 2026-11-19, `fix/panaroo-suffixes`]** `resolve_locus()` (same code in all three scripts) tries the raw ID, then strips `_len`/`_pseudo`/`_stop` one at a time (the user's real Panaroo run has `_stop`, `_len` and `_len_stop`). Test data vs baseline v3: missing loci 130 → 0; 41/42/47 more genes per isolate; 53 new rows; core 5674 → 5713; no existing kallisto value changed; 6 raw/norm values in 4 rows moved (median length gained an isolate); size factors unchanged; nTPM × 0.997163 uniformly (global constant only); UMAP identical. Baseline re-captured as **v4**. 4 new unit tests (they fail on the old code); 40/40 pass in the image.
- **R8. Minor issues** [fixed 2026-11-18: unused `gpa` input, `ch_trimgalore_*` channels, `when: metadata.exists()` and `custom_runName` removed; `--group` moved to optional args in help; undefined `params.workdir` removed from the slurm profile; conda profile now sets `conda.enabled = true` (the old `params.enable_conda` did nothing)] (fine to fix with no output change): `KALLISTO_QUANT` takes the `gpa`
  input but never uses it; `ch_trimgalore_*` channels are never consumed; `when: metadata.exists()`
  can never be false; the `--group` help text says it's mandatory and also gives a default; the
  `slurm` profile needs `params.workdir`, which is never defined.

## Open questions for the user

- **Q15 (open, 2026-11-19):** merge the DOI docs commit (078b454, docs only) from `dev` into
  `main` now, so the GitHub front page shows the DOI, or leave it for the next release?
- Q1–Q3: resolved 2026-11-18 (see Decision log).
- Q1b: resolved 2026-11-18. Keep `ntpm.tsv`; the "n" (normalised) is explained in the README and docs.
- Q9: resolved 2026-11-18 → P9 (fixed seed; rewrite the UMAP in Python).
- Q10: resolved 2026-11-18 → P10 (`--strandedness auto`).
- Q11: resolved 2026-11-18. `--strandedness auto` is the default; the test profile uses it.
- Q12: resolved 2026-11-19. Yes: strip Panaroo status suffixes when matching locus IDs (done, R12).
- Q5: resolved 2026-11-18. Name is **PanScript**.
- Q14: resolved 2026-11-19 (user: go ahead). Done in H15. Was: raise the minimum Nextflow from 23.10 to 24.10? It enables `resourceLimits` (replacing the `max_cpus`/`max_memory` params and closures) and a current nf-schema (2.4.x; 2.6+ needs 25.10 and Java 17). Also: right-size the process labels (kallisto currently asks for 12 CPUs / 72 GB, far more than bacterial gene sets need).
- Q13: resolved 2026-11-19. Adopt all the cheap wins in `docs/nf-core-assessment.md` (H9–H14); the user approved the output changes they bring.
- Q4: resolved 2026-11-19. Reject mixed single-end and paired-end input with a clear error (R5).
- Q6: resolved 2026-11-18 → P7. The FASTA files are per-isolate CDS from Bakta or Prokka, the same
  annotation that goes into Panaroo, supplied by the user. Pseudogenes are wanted *if* the user
  supplies their sequences.
- Q7b resolved → P8. Q7c: the user should check any real analyses of this cohort for the
  strandedness setting they used.
- Q8: resolved 2026-11-18. Fail when more than 50% of a sample's loci are missing from its
  FASTA; warn otherwise.
- Q7a: resolved 2026-11-18. The test profile now uses `strandedness = 1`, and the baseline was
  re-captured.

## Decision log

| Date | Decision | Reasoning |
|------|----------|-----------|
| 2026-11-18 | Created CLAUDE.md and PROGRESS.md for multi-agent coordination. | User request. |
| 2026-11-18 | Remove `--perc`; output all genes with NA for absent genes. | It was meant as a gene-inclusion filter but was never wired in; not needed. |
| 2026-11-18 | Drop TMM; DESeq2 median-of-ratios only. | Simpler; avoids the TMM scale inconsistency (R3). |
| 2026-11-18 | The length-normalised measure must use core-gene scaling, not per-sample totals over all genes. | Different accessory content between isolates would distort per-sample totals. |
| 2026-11-18 | Keep split-gene averaging of the per-isolate rate (pending Q2 on median length). | Mean/mean = sum/sum, so the rate is correct. |
| 2026-11-18 | Length-normalised measure = option (b): (count / own effective length) / DESeq2 core size factor × global constant to TPM-like scale. Rejected plain TPM and core-anchored TPM. | Robust to differing accessory content and to a few dominant core genes. Comparable between genes and across samples. |
| 2026-11-18 | Split genes: sum counts and lengths (P5). | Median length then reflects the intact gene. |
| 2026-11-18 | Output linear matrices, not log2 (P6). | More conventional; users log-transform as needed. Must be stated clearly. |
| 2026-11-18 | Loci missing from user FASTA: quantify what's provided, keep missing ones NA, and report them in a log (P7). | Don't silently discard genes, but sequences are the user's responsibility. |
| 2026-11-18 | Add an automatic strandedness check (P8). | The test data runs with the wrong strandedness (R10); this protects users from the same mistake. |
| 2026-11-18 | Phase 1 implemented on `feat/normalisation`: P4, P3, P6, P5, P2, each checked against the baseline and unit tests. | See task rows for the observed differences. |
| 2026-11-18 | Keep the name `ntpm.tsv` (nTPM). | Confirmed by user; explain the "n" in the README and docs. |
| 2026-11-18 | UMAP: fixed seed, rewritten in Python (P9). | Reproducible figures; the other scripts are mostly Python. |
| 2026-11-18 | Order: Phase 1b (P7, P8) next on the current container, then Phase 2 bundles P1 (fastp) + P9 (Python UMAP) + H2 (pinned rebuild) into one container rebuild. | P7/P8 need no new tools, so they can be checked cleanly against baseline v2. Every tool-version change is then confined to one rebuild and one baseline re-capture. |
| 2026-11-18 | Phase 2 (fastp + Python UMAP + pinned container) stays on `feat/phase2-container` until the new image is published. | Merging earlier would point `dev` at the old `:latest` image, which lacks fastp and umap-learn, so docker runs of `dev` would break. |
| 2026-11-18 | Published `adamd3/strainseq:2.0.0-dev`; `--strandedness auto` is the default. | User approved both. The name stays as is for now (H1 open); the image can be re-tagged at rename. |
| 2026-11-18 | Pipeline renamed to **PanScript**; image `adamd3/panscript`. | User decision after several rounds (see Q5). Docker needs lowercase names. |
| 2026-11-18 | Git workflow: `main` (releases) / `dev` (integration) / feature branches; no AI attribution in commits. | User delegated branching; user request on attribution. |
| 2026-11-19 | Strip Panaroo status suffixes when matching locus IDs (Q12/R12); baseline v4. | The suffixes aren't part of the sequence ID, so ~0.6% of genes per isolate (39 core) were NA for no biological reason. Raw ID is tried first so a real ID ending in `_stop` still wins. |
| 2026-11-19 | Released v2.0.0 (tag, GitHub release, image `adamd3/panscript:2.0.0`, Zenodo DOI 10.5281/zenodo.24375403). | User approved the image push and release. |
| 2026-11-19 | Validate params with our own schema-driven code, not the nf-schema plugin (H9). | nf-schema 2.4.2 breaks on Nextflow 26 (string CLI values, ignored config); unpinned versions print a warning banner on every run. Own code is version-proof and needs no plugin download. |
| 2026-11-19 | Minimum Nextflow raised to 24.10 (H15). | Needed for `resourceLimits` and a current nf-schema; 24.10 is two years old. Users on 23.x must upgrade. |
| 2026-11-19 | No site-specific profiles in `nextflow.config` (H14). | User: the `cambridge`/`slurm` profiles are relics; site settings belong in a user `-c` config. |
| 2026-11-19 | Reject mixed single-end/paired-end runs (Q4/R5). | Simplest correct behaviour; kallisto runs in one mode for all samples. Users split such cohorts into separate runs. |

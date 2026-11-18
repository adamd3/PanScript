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

## Plan (overview)

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

**Phase 2: Trimming (P1, fastp), on `feat/fastp`.** This can start in parallel with Phase 1. It
needs fastp in the container: rebuild the image or use a biocontainer per process. Expected:
counts change slightly (different trimmer), so compare correlations with the baseline rather than
requiring identical values.

**Phase 3: Reproducibility and repo hygiene.** H2 (pin versions and version-tag the image), H3
(fix CI, plus a test option that reads local files), H4/H7/H8 (cleanup), H5 (README), then
H6 (nf-core cost/benefit write-up).

**Phase 4: Name and release (H1).** Rename everything. Merge `dev` into `main` and tag `v1.0.0`.
Get a Zenodo DOI.

**Phase 5: Manuscript (M1–M3).** Outline and benchmarking design can start any time, but the
results need the Phase 1–2 pipeline.

---

## Tasks

### Pipeline changes (these affect outputs, so get user sign-off on the details)

| # | Task | Status | Notes |
|---|------|--------|-------|
| P1 | Replace Trim Galore with fastp for adaptor trimming; drop FastQC (fastp does QC). | TODO | Approved. Also update the `KALLISTO_QUANT` filename symlinks (they rely on Trim Galore's names), `environment.yml`, README and `save_trimmed` publishing. Consider publishing the fastp JSON/HTML reports, and optionally MultiQC (already in env). Agree fastp parameters with the user. |
| T0 | Capture baseline outputs and add a regression check. | DONE (2026-11-18, `feat/regression-baseline`, merged into `dev`) | The baseline is in `tests/baseline/` (provenance in `PROVENANCE.md`), and `tests/compare_to_baseline.py` checks new runs against it (it matches the baseline itself exactly). Confirmed R1, R2, R4 and the log2 outputs on real data. Found **R9** and **R10**. Still to do: synthetic-data unit tests (T1). |
| T1 | Unit tests on synthetic data for the merge and normalisation scripts. | IN PROGRESS (claude, 2026-11-18) | `tests/unit/test_normalisation.py`: 7 tests, all passing against *current* behaviour, using an independent numpy reference. Update the expected values with each Phase 1 step. | The test dataset is too sparse to check normalisation maths (see R10). Build tiny hand-computable inputs (split genes, accessory NA, known size factors) and check each output formula, especially nTPM. Write it before or during Phase 1. |
| P2 | Replace RPKM with nTPM, computed from DESeq2 core-gene size factors. | TODO | Approved 2026-11-18 (option b). Definition: x_gs = (c_gs / L_gs) / sf_s, where L is the isolate's own effective length, then multiplied by one global constant K so that the mean over samples of the per-sample total (all non-NA genes) is 1e6. K is identical for every gene and sample, so it doesn't affect comparisons. Proposed name: **nTPM** ("normalised TPM", a precedent from the Human Protein Atlas), file `ntpm.tsv`; check the name with the user. The README and methods must state that per-sample totals are only *approximately* 1e6. Median length cancels out, so P5 doesn't affect this output. |
| P3 | Remove `--perc`, `SUBSET_GENES` and `gene_set_ST_*.tsv`; output every gene. | TODO | Approved 2026-11-18. Also remove the `-p` branch in the R scripts, `bin/subset_genes.py`, and the param from config, help, README and validation. Accessory genes stay NA. |
| P4 | Remove TMM: DESeq2 only. | TODO | Approved 2026-11-18. Remove `TMM_NORMALISE_COUNTS`, `bin/TMM_normalise_counts.R` and `--norm_method` (help, config, README, validation). edgeR can then come out of `environment.yml`. |
| P5 | Sum (not average) counts and effective lengths for split genes. | TODO | Approved 2026-11-18. Edit `bin/merge_kallisto_counts.py` and `bin/merge_kallisto_lens.py`. The per-isolate rate is unchanged; the median length, and so the magnitude of `raw_counts` and `norm_counts`, changes only for genes that are split in some isolates. |
| P6 | Output linear (non-log) matrices. | TODO | Approved 2026-11-18. Drop the hard-coded `-t TRUE`, and say clearly in the README (and in file names or headers if useful) that values are untransformed. **`UMAP_SAMPLES` currently relies on log2 `norm_counts`, so `umap.R` must apply log2(x+1) itself** to keep the plot unchanged. |
| P7 | Report loci that Panaroo lists but that are missing from the isolate FASTA. | TODO | Approved 2026-11-18 (R9/Q6): use whatever sequences the user provides, and **report** what's missing. Don't change values (missing loci stay NA, as now). Proposed design: a per-sample check before kallisto compares the Panaroo loci for that `dna_sample_id` with the FASTA headers and writes `pipeline_info/missing_loci.tsv` (sample, gene group, locus, flags such as `refound_stop` and `partial_split`, meaning only some loci of a split cell are present), plus a warning with counts per sample. Proposal (needs user OK, since it can stop runs): **fail** if most loci are missing (for example > 50%), which almost certainly means the wrong FASTA or a mismatched `dna_sample_id`. Note that with P5 (summing), a `partial_split` gene undercounts, so it must be flagged. |
| P8 | Automatic strandedness check. | TODO | Approved 2026-11-18 (Q7b). Quantify a subsample of reads (for example 200k) per sample in both orientations (kallisto `--fr-stranded` / `--rf-stranded`). Infer forward/reverse/unstranded from the ratio, and **warn** when it disagrees with `--strandedness` (and possibly offer `--strandedness auto` to use the inferred value, as nf-core/rnaseq does). Report the result per sample in `pipeline_info/`. Thresholds and the auto mode need user sign-off. |
| — | Suggested order: T0 → one branch `feat/normalisation` for P3 + P4 + P5 + P6 + P2 (they all touch the same R and Python scripts) → P1 (fastp; independent, can run in parallel on its own branch) → H-tasks. | | |

### Publication / repo hygiene (no change to outputs)

| # | Task | Status | Notes |
|---|------|--------|-------|
| H1 | Choose a new pipeline name. | NEEDS USER (pick from candidates) | History: OrthoQuant (original) → StrainQuant. The GitHub repo is currently `adamd3/StrainQuant`, and the old names `OrthoQuant` and `StrainSeq` redirect to it. Renaming the repo back to OrthoQuant would work smoothly with GitHub's redirects. Candidates proposed 2026-11-18: **OrthoQuant** (recommended), **PanQuant**. A quick web search found no tool with either exact name; the nearest are *OrthoQuantum* (eukaryotic phylogenetic profiling, NAR 2022) and *PanXpress* (bacterial pan-transcriptomics). Before deciding, check bioconda, PyPI and GitHub for clashes. It should describe what the pipeline does: per-isolate pseudo-alignment and pangenome-aware merging of expression. It shouldn't imply "strains". Then rename everything: repo, README, `main.nf` header and help, Docker image `adamd3/strainseq`, conda env name, CI, test URLs. |
| H2 | Pin software versions and rebuild the container. | TODO | **Not before T0; combine the rebuild with P1** (fastp has to go in the image anyway). Dockerfile/env issues: (1) base `continuumio/miniconda3` is unpinned; pin it or switch to a pinned `mambaorg/micromamba`. (2) `environment.yml` uses `>=`; pin exact versions. (3) Drop the Anaconda `defaults` channel (licensing terms for large organisations); use conda-forge + bioconda with strict priority. (4) `COPY . /app` bakes the whole repo (test data, stray files) into the image; it isn't needed, since Nextflow supplies `bin/`. Remove it or add a `.dockerignore`. (5) Update dependencies: add fastp; remove trim-galore, fastqc and (after P4) edgeR; keep multiqc only if used. (6) Replace `:latest` and `--pull=always` with a version tag `adamd3/<name>:vX.Y.Z`, set once in config (modules also hard-code the image). (7) Build for amd64 **and** arm64 with `docker buildx` *if* all bioconda dependencies have linux-aarch64 builds; check each one, otherwise stay amd64-only. After rebuilding, compare outputs with the T0 baseline: a new kallisto or DESeq2 version could shift values slightly. |
| H3 | Fix CI. | TODO | No run history shows on GitHub. That's expected: Actions is enabled and the secrets exist, but the last push was in July 2025, so old runs have passed the retention limit. The push job's `if:` checks `github.repository == 'adamd3/StrainQuant'`. That's correct today, but it must be updated if the repo is renamed (H1). **Side effect to fix before using PRs:** every PR run rebuilds the image from unpinned `environment.yml` and **pushes it as `adamd3/strainseq:latest`**. That silently changes the tools that every `-profile docker` run (and the T0 baseline digest) uses. Until H3 is fixed, merge branches into `dev` locally instead of through GitHub PRs. It pushes `:latest` on every PR, and fork PRs have no secrets. It only triggers on pushes to `main` (add `dev`). It uses old actions (`checkout@v2`, `login-action@v1`). Test data comes from `main`, not the branch under test. The minimum Nextflow version is 21.10.3, which is old. Consider adding small unit tests (pytest/testthat) for the merging and normalisation scripts. |
| H4 | Clean up unused code. | TODO | See "Unused files" in CLAUDE.md. Confirm with the user before deleting (some may be useful as helper utilities). |
| H5 | Update README. | TODO | Use the agreed definitions of the output matrices (CLAUDE.md, "Meaning of the output matrices"). Update parameters after P3/P4. Add a citation section. Do this after P1–P6 land. |
| H6 | Assess nf-core compliance. | TODO | The user is reluctant to do heavy lifting just for nf-core. First write up a cost/benefit: what nf-core listing needs (template, nf-schema, modules, MultiQC, nf-test) versus cheap wins worth taking anyway (`nextflow_schema.json`, `-profile test`, versioned containers, `versions.yml`). |
| H7 | Clean up the repo. | TODO | Stray zero-byte `test_data/.!NNNNN!*.fq.gz` files (failed-transfer leftovers, untracked). Decide whether to commit `docker_build_and_push.sh`. Add a CITATION.cff and a Zenodo DOI before submission. |
| H8 | Fix minor code issues (R8). | TODO | No output change. |

### Manuscript (Genome Biology)

| # | Task | Status | Notes |
|---|------|--------|-------|
| M1 | Draft an outline for a Genome Biology "Method" or "Software" article. | TODO | Check the current GB article types and requirements: availability, OSI licence (MIT ✓), archived version. |
| M2 | Plan benchmarking and validation. | TODO | **Related tools to compare against or discuss:** PanXpress (PyPI `panxpress`: pan-transcriptomic gapped k-mer index, alignment-free, includes its own annotation-correction step); the Greene lab "2022-microberna" single-pangenome versus best-genome comparison (github.com/greenelab/2022-microberna); single-reference mapping. Run a full literature search during M1. For example, compare against single-reference mapping (such as PAO1) on multi-isolate datasets; use simulated reads with known expression for different isolate sets and accessory fractions; check robustness of core-gene size factors as core size shrinks. The design needs input from the user. |
| M3 | Choose the showcase dataset(s). | TODO | Test data is *P. aeruginosa* (SRX5123741/2/4, the HZI AMR cohort). |

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
  mixed input be supported, or rejected with a clear error?
- **R6. Log2 outputs.** `norm_counts.tsv` and `rpkm_counts.tsv` are log2(x+1)-transformed
  (`-t TRUE` is hard-coded), but the file names and README don't say so (P6).
- **R7. Downstream differential expression.** Length scaling gives non-integer `raw_counts`. If
  users take them into DESeq2 for differential expression, document the recommended route:
  rounded `raw_counts`, or raw `est_counts` plus a per-gene, per-sample length offset (as tximport
  does). Unanswered; lower priority, but reviewers may ask.
- **R9. Pseudogenised "refound" genes become NA.** Found 2026-11-18 in the baseline. Panaroo
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
- **R8. Minor issues** (fine to fix with no output change): `KALLISTO_QUANT` takes the `gpa`
  input but never uses it; `ch_trimgalore_*` channels are never consumed; `when: metadata.exists()`
  can never be false; the `--group` help text says it's mandatory and also gives a default; the
  `slurm` profile needs `params.workdir`, which is never defined.

## Open questions for the user

- Q1–Q3: resolved 2026-11-18 (see Decision log).
- **Q1b (P2):** Is the name "nTPM" / `ntpm.tsv` OK?
- **Q4 (R5):** Mixed single-end and paired-end input: support it or reject it?
- **Q5 (H1):** Name ideas: should the name stress pangenome, isolate-specific or reference-free?
- Q6: resolved 2026-11-18 → P7. The FASTA files are per-isolate CDS from Bakta or Prokka, the same
  annotation that goes into Panaroo, supplied by the user. Pseudogenes are wanted *if* the user
  supplies their sequences.
- **Q7a (R10):** Set `strandedness = 1` in the test profile? That changes the baseline, so it
  would be re-captured. (Q7b resolved → P8. Q7c: the user should check any real analyses of this
  cohort for the setting they used.)
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
| 2026-11-18 | Git workflow: `main` (releases) / `dev` (integration) / feature branches; no AI attribution in commits. | User delegated branching; user request on attribution. |

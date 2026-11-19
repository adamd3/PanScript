# nf-core: cost/benefit (PROGRESS.md H6)

Written 2026-11-19 to help decide whether PanScript should aim for an nf-core listing. nf-core's
requirements change with each template release, so check <https://nf-co.re/docs/guidelines> before
acting on the details below.

## Recommendation

**Don't aim for an nf-core listing before the Genome Biology submission.** Instead, take the
cheap wins listed below, which give most of the user-facing benefit (parameter validation, a
run report, recorded tool versions, a familiar layout) for a few days' work. Keep the code close
enough to nf-core conventions that a port stays feasible later, if the user base justifies it.

Genome Biology doesn't require nf-core. Reviewers will look for what PanScript already has or is
about to get: versioned containers, CI, tests, a tagged release and a Zenodo DOI.

## What a listing would require

| Requirement | Where PanScript is now | Effort |
|---|---|---|
| Proposal accepted by the nf-core community (Slack `#new-pipelines`), with no overlap with an existing pipeline. The repo moves to the `nf-core` GitHub organisation as `nf-core/panscript`. | Own repo. Overlap is probably low (nf-core/rnaseq is reference-based), but that's for the community to decide. | Low, but outside our control. Ownership is shared with nf-core afterwards. |
| Rebuild on the nf-core template (`nf-core pipelines create`) and keep it in sync with each template release (sync PRs several times a year). | Custom layout (`main.nf`, `modules/`, `configs/`). | **High** to start, then ongoing. |
| One container per process from Bioconda/Biocontainers (or Seqera Containers), not a single image. The `bin/` scripts need containers with pandas, R/DESeq2 and umap-learn. | One pinned image (`adamd3/panscript`). | **Medium–high.** Means re-validating against the baseline. |
| Use nf-core/modules where they exist (`fastp`, `kallisto/index`, `kallisto/quant`); local modules follow the module spec (`meta` maps, `ext.args`, `versions.yml` or topic channels). | Local modules with a different interface. | Medium. |
| `nextflow_schema.json` + the nf-schema plugin; CSV samplesheet validated against `assets/schema_input.json`. | TSV metadata with free-form extra columns (used for `--group`); validation in `check_metadata.py`. | Medium. The samplesheet format changes for users. |
| A MultiQC report. | None (fastp JSON reports only). MultiQC isn't in the image. | Low–medium. |
| nf-test pipeline and module tests with snapshots; `nf-core pipelines lint` passing in CI; `test` and `test_full` profiles (full-size test on AWS). | Python unit tests + baseline comparison in CI; `test` profile. | Medium. Duplicates what our tests already cover. |
| Docs in nf-core format (`docs/usage.md`, `docs/output.md`, `CITATIONS.md`, `CHANGELOG.md`); site configs go in nf-core/configs. | README + `docs/normalisation.md`; a `cambridge` profile in `nextflow.config`. | Low. |
| MIT licence. | MIT. | None. |

Overall, roughly **several weeks** of work, plus a full re-validation, and ongoing maintenance
through template syncs and reviews.

## What a listing would give

- Visibility: listed on nf-co.re, discoverable by the nf-core community.
- Credibility, and a shared support channel (nf-core Slack).
- Institutional profiles for free (`-profile <institution>` via nf-core/configs).
- Launchable from Seqera Platform with a generated form (also available from the schema alone).

Against that: shared ownership, the template-sync burden, and a citation split. nf-core pipelines
are typically cited via the nf-core paper plus the pipeline DOI, alongside our own paper.

## Cheap wins worth taking anyway

In rough order of value for effort. Items marked **(output)** change pipeline outputs, so they need
the user's OK (CLAUDE.md rule 1).

1. **`nextflow_schema.json`** for parameter validation and generated `--help`, and so the
   pipeline is launchable from Seqera Platform. *Done 2026-11-19 (H9), validated by our own code in
   `main.nf` rather than the nf-schema plugin, which misbehaves on Nextflow ≥ 25 at the version
   that supports our 24.10 minimum.*
2. **Record tool versions** in `pipeline_info/software_versions.yml` (kallisto, fastp, R, DESeq2,
   Python, pandas, umap-learn). This is useful for any methods section. **(output: new file)**
3. **`--publish_dir_mode`, default `copy`.** Most outputs are currently symlinks into `work/`, so
   deleting `work/` silently breaks `results/` (a common user trap). **(output: files become
   copies)**
4. **MultiQC report** summarising fastp and kallisto (MultiQC supports both). It needs MultiQC in
   `environment.yml`, so do it with the next image rebuild (for example the v2.0.0 release
   image). **(output: new report)**
5. **nf-core-style docs layout:** `docs/usage.md`, `docs/output.md`, `CITATIONS.md`,
   `CHANGELOG.md`. Most of this is already planned for the v2.0.0 release.
6. **Move the `cambridge` profile** to a separate config file (or `-c` at run time), so
   `nextflow.config` holds only general profiles.

Not worth doing outside a listing: per-process Biocontainers, nf-core/modules, nf-test (our unit
tests + baseline comparison already cover this) and the full template.

## If the decision is revisited later

Keeping process interfaces simple (one tool or script per process, `meta`-style tuples where
cheap) and adopting items 1–5 above would make a later port mostly mechanical. The main remaining
work would be the template, per-process containers and the samplesheet format change.

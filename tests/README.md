# Tests

## Unit tests (synthetic data)

```bash
python3 -m unittest discover -s tests/unit -v
# or, with exactly the pipeline's software versions (needed for the UMAP tests):
docker run --rm -v $PWD:/repo -w /repo <image> python -m unittest discover -s tests/unit
```

`tests/unit/test_normalisation.py` builds a tiny hand-made dataset (core genes with known
library-depth differences, accessory genes, a split gene). It runs the `bin/` merge and
normalisation scripts on it and compares every output with an independent numpy implementation of
the formulas in `docs/normalisation.md`. It needs python3 with pandas, and Rscript with DESeq2,
optparse, matrixStats and tidyverse (all present in the pipeline container). **When a
normalisation formula changes, update the reference implementation in the test as part of the same
change.**

# Regression tests

`baseline/` holds the outputs of the `-profile test` run on unchanged code, as a fixed record of
pipeline behaviour before the 2026 publication changes (PROGRESS.md task T0). `baseline/PROVENANCE.md`
records exactly how they were produced: commit, Nextflow version and container digest.

`baseline_pe/` is the same for `-profile test_pe`: paired-end reads simulated by
`simulate_pe_reads.py` (deterministic; re-running it reproduces `test_data/pe/` exactly) with a
known truth. Besides the baseline comparison (`--baseline tests/baseline_pe`),
`check_pe_truth.py <results>` checks the inferred strandedness and that counts track the
simulated truth. CI runs both test profiles.

## Checking a change against the baseline

```bash
nextflow run . -profile test,docker --outdir <results> -w <work>
python3 tests/compare_to_baseline.py <results>
```

The script compares each matrix cell by cell (NA positions must match; numbers must agree within
`--rtol`/`--atol`) and lists the rows that changed most. If a change is *meant* to alter an output
in a known way, transform the baseline so the comparison is still exact. For example, after
outputs switch from log2(x+1) to linear:

```bash
python3 tests/compare_to_baseline.py <results> \
    --transform norm_counts=exp2m1 --only raw_counts norm_counts
```

Each Phase 1 step in PROGRESS.md states which differences it should produce. Anything else is a
bug, or something to raise with the user.

The baseline is replaced only deliberately, for example after the container is rebuilt (H2). Keep
the old one, or record why it was replaced.

# Paired-end baseline provenance

Reference outputs for the paired-end test profile (`test_pe`), checked in CI with
`tests/compare_to_baseline.py <results> --baseline tests/baseline_pe` and, against the
simulation truth, `tests/check_pe_truth.py <results>`.

| Item | Value |
|------|-------|
| Date | 2026-11-20 (PE baseline v1) |
| Pipeline code | `feat/pe-test-data` (branched from `dev` @ `46a4cea`) |
| Command | `nextflow run . -profile test_pe,docker --outdir <results>` (no overrides) |
| Nextflow | 26.10.6 (macOS arm64; container under amd64 emulation) |
| Container | `adamd3/panscript:2.0.0` (`adamd3/panscript@sha256:7a1fbe6c1d3dd916f37d60b515fb04c9253e785b0f38f5646e7e3443ff49422b`) |
| Inputs | `test_data/pe/`, simulated by `tests/simulate_pe_reads.py` (deterministic): 4 samples (PE_A, PE_B, PE_C1, PE_C2; C1/C2 are replicates of isolate SRR8737283), 30–60k read pairs, 2 x 100 bp, reverse-stranded, fragments ~ N(200, 50) |
| Parameters | defaults (`strandedness=auto`, fastp defaults), `group=isolate` |

## Known properties of this baseline
- Strandedness inferred reverse (2) for all samples (rf 99.7–99.99%, fr 0.03–0.27%).
- kallisto: ≥ 99.98% of pairs pseudo-aligned (simulated reads, no contamination); estimated
  mean fragment length 198.4–199.4 (simulated 200).
- fastp trimmed adapters from about 2% of reads (fragments shorter than 100 bp); no reads removed.
- 3 kallisto indexes for 4 samples (one per isolate).
- Versus truth: Pearson r (log1p, loci or groups with ≥ 5 fragments) 0.997–0.999.

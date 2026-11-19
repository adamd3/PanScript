#!/usr/bin/env python3
"""Unit tests for bin/check_metadata.py (paired-end detection, PROGRESS.md R5/Q4)."""

import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "check_metadata.py"


class CheckMetadataTest(unittest.TestCase):
    def setUp(self):
        self.d = Path(tempfile.mkdtemp(prefix="sq_meta_"))

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def run_script(self, fastq2):
        n = len(fastq2)
        pd.DataFrame({"rna_sample_id": [f"R{i}" for i in range(n)],
                      "dna_sample_id": [f"D{i}" for i in range(n)],
                      "sample_name": [f"S{i}" for i in range(n)],
                      "fastq1": [f"r{i}_1.fq.gz" for i in range(n)],
                      "fastq2": fastq2,
                      "fasta": [f"d{i}.fna" for i in range(n)]}).to_csv(
            self.d / "meta.tsv", sep="\t", index=False)
        return subprocess.run([sys.executable, str(SCRIPT), "meta.tsv", "out.tsv"],
                              cwd=self.d, capture_output=True, text=True)

    def paired_column(self):
        return pd.read_csv(self.d / "out.tsv", sep="\t", dtype=str)["paired"].tolist()

    def test_all_single_end(self):
        r = self.run_script(["", "", ""])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.paired_column(), ["0", "0", "0"])

    def test_all_paired_end(self):
        r = self.run_script(["r0_2.fq.gz", "r1_2.fq.gz"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.paired_column(), ["1", "1"])

    def test_mixed_is_rejected(self):
        r = self.run_script(["r0_2.fq.gz", "", "r2_2.fq.gz"])
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("mixes single-end samples (no fastq2: S1)", r.stderr)
        self.assertIn("paired-end samples (fastq2 given: S0, S2)", r.stderr)


if __name__ == "__main__":
    unittest.main()

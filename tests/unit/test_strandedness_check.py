#!/usr/bin/env python3
"""Unit tests for bin/strandedness_check.py (PROGRESS.md P8)."""

import gzip
import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "strandedness_check.py"
spec = importlib.util.spec_from_file_location("strandedness_check", SCRIPT)
sc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc)


class ClassifyTest(unittest.TestCase):
    def test_thresholds(self):
        self.assertEqual(sc.classify(900, 100)[0], 1)    # f = 0.9 forward
        self.assertEqual(sc.classify(800, 200)[0], 1)    # boundary 0.8
        self.assertEqual(sc.classify(100, 900)[0], 2)    # reverse
        self.assertEqual(sc.classify(200, 800)[0], 2)    # boundary 0.2
        self.assertEqual(sc.classify(500, 500)[0], 0)    # unstranded
        self.assertIsNone(sc.classify(700, 300)[0])      # ambiguous
        self.assertIsNone(sc.classify(0, 0)[0])          # nothing aligned

    def test_test_dataset_values(self):
        # observed on the test data (PROGRESS.md R10): forward-stranded
        self.assertEqual(sc.classify(50508, 6059)[0], 1)
        self.assertEqual(sc.classify(50003, 2817)[0], 1)


class ScriptTest(unittest.TestCase):
    def setUp(self):
        self.d = Path(tempfile.mkdtemp(prefix="sq_strand_"))

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def run_script(self, *args):
        return subprocess.run([sys.executable, str(SCRIPT), *args], cwd=self.d,
                              capture_output=True, text=True)

    def write_fq(self, name, n):
        with gzip.open(self.d / name, "wt") as fh:
            for i in range(n):
                fh.write(f"@r{i}\nACGT\n+\nIIII\n")

    def test_subsample_paired_and_short_input(self):
        self.write_fq("a_1.fq.gz", 10)
        self.write_fq("a_2.fq.gz", 10)
        r = self.run_script("subsample", "--n", "4", "--out", "sub", "a_1.fq.gz", "a_2.fq.gz")
        self.assertEqual(r.returncode, 0, r.stderr)
        for i in (1, 2):
            lines = gzip.open(self.d / f"sub_{i}.fq.gz", "rt").read().splitlines()
            self.assertEqual(len(lines), 16)
            self.assertEqual(lines[0], "@r0")
        # fewer reads than requested: keeps them all
        r = self.run_script("subsample", "--n", "100", "--out", "s2", "a_1.fq.gz")
        self.assertEqual(len(gzip.open(self.d / "s2_1.fq.gz", "rt").read().splitlines()), 40)

    def infer(self, fr, rf, specified):
        for mode, n in [("fr", fr), ("rf", rf)]:
            (self.d / f"{mode}.json").write_text(json.dumps(
                {"n_processed": 1000, "n_pseudoaligned": n}))
        r = self.run_script("infer", "--fr", "fr.json", "--rf", "rf.json", "--sample", "S1",
                     "--specified", str(specified), "--out", "out.tsv")
        self.assertEqual(r.returncode, 0, r.stderr)
        return r, pd.read_csv(self.d / "out.tsv", sep="\t", keep_default_na=False).iloc[0]

    def test_agreement_is_silent(self):
        r, row = self.infer(500, 50, specified=1)
        self.assertEqual(row["inferred"], "forward")
        self.assertEqual(row["agrees"], "yes")
        self.assertEqual(r.stdout.strip(), "")

    def test_mismatch_warns_with_suggestion(self):
        r, row = self.infer(500, 50, specified=2)
        self.assertEqual(row["agrees"], "NO")
        self.assertIn("--strandedness 1", r.stdout)

    def test_used_strandedness_is_specified_value(self):
        _, row = self.infer(500, 50, specified=2)       # mismatch still uses 2
        self.assertEqual(row["used_strandedness"], 2)

    def test_auto_uses_inferred_value(self):
        r, row = self.infer(50, 500, specified="auto")
        self.assertEqual(row["inferred"], "reverse")
        self.assertEqual(row["used_strandedness"], 2)
        self.assertEqual(row["agrees"], "n/a (auto)")
        self.assertEqual(r.stdout.strip(), "")

    def test_auto_undetermined_is_an_error(self):
        for mode, n in [("fr", 700), ("rf", 300)]:
            (self.d / f"{mode}.json").write_text(json.dumps(
                {"n_processed": 1000, "n_pseudoaligned": n}))
        r = self.run_script("infer", "--fr", "fr.json", "--rf", "rf.json", "--sample", "S1",
                            "--specified", "auto", "--out", "out.tsv")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("could not determine", r.stderr)
        # the table is still written for inspection
        self.assertTrue((self.d / "out.tsv").exists())

    def test_undetermined_warns(self):
        r, row = self.infer(700, 300, specified=1)
        self.assertEqual(row["inferred"], "undetermined")
        self.assertEqual(row["agrees"], "unknown")
        self.assertIn("could not determine", r.stdout)


if __name__ == "__main__":
    unittest.main()

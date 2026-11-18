#!/usr/bin/env python3
"""Unit tests for bin/check_fasta_loci.py (PROGRESS.md P7)."""

import gzip
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "check_fasta_loci.py"


def write_fasta(path, ids, gz=False):
    text = "".join(f">{i} some annotation\nACGT\n" for i in ids)
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)


class CheckFastaLociTest(unittest.TestCase):
    def setUp(self):
        self.d = Path(tempfile.mkdtemp(prefix="sq_loci_"))
        # isolate dA: 10 gene groups
        cells = {f"g{i}": f"dA_{i}" for i in range(1, 8)}
        cells["refound"] = "1_refound_99_stop"
        cells["refoundL"] = "1_refound_77_len"
        cells["splitP"] = "dA_sp1;dA_sp2"           # one locus missing -> partial
        cells["splitA"] = "dA_sa1;dA_sa2"           # both missing
        cells["absent"] = ""                        # not in this isolate
        rows = [{"Gene": g, "Non-unique Gene name": "", "Annotation": "", "dA": c, "dB": "x"}
                for g, c in cells.items()]
        pd.DataFrame(rows).to_csv(self.d / "gpa.csv", index=False)
        self.present = [f"dA_{i}" for i in range(1, 8)] + ["dA_sp1"]

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def run_script(self, fasta, dna="dA", max_frac=0.5):
        r = subprocess.run(
            [sys.executable, str(SCRIPT), "--gene_presence_absence", "gpa.csv",
             "--dna_sample_id", dna, "--sample_name", "S1", "--fasta", fasta,
             "--max_missing_frac", str(max_frac), "--outf", "miss.tsv",
             "--summary", "summary.tsv"],
            cwd=self.d, capture_output=True, text=True)
        return r

    def test_reports_missing_loci_with_flags(self):
        write_fasta(self.d / "a.fna", self.present)
        r = self.run_script("a.fna")
        self.assertEqual(r.returncode, 0, r.stderr)
        miss = pd.read_csv(self.d / "miss.tsv", sep="\t").fillna("")
        got = dict(zip(miss["locus"], miss["flags"]))
        self.assertEqual(got, {"1_refound_99_stop": "refound_stop",
                               "1_refound_77_len": "refound_len",
                               "dA_sp2": "split_partial",
                               "dA_sa1": "split_all_missing",
                               "dA_sa2": "split_all_missing"})
        s = pd.read_csv(self.d / "summary.tsv", sep="\t").iloc[0]
        self.assertEqual(s.loci_listed, 13)
        self.assertEqual(s.loci_missing_from_fasta, 5)
        self.assertEqual(s.gene_groups_na_due_to_missing, 3)   # refound, refoundL, splitA
        self.assertEqual(s.gene_groups_split_partial, 1)
        self.assertEqual(s.loci_refound_stop, 1)
        self.assertEqual(s.loci_refound_len, 1)
        self.assertIn("5 of 13 loci", r.stdout)

    def test_no_message_when_nothing_missing(self):
        write_fasta(self.d / "a.fna", self.present + ["1_refound_99_stop", "1_refound_77_len", "dA_sp2",
                                                       "dA_sa1", "dA_sa2"])
        r = self.run_script("a.fna")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(r.stdout.strip(), "")
        self.assertEqual(len(pd.read_csv(self.d / "miss.tsv", sep="\t")), 0)

    def test_gzipped_fasta(self):
        write_fasta(self.d / "a.fna.gz", self.present, gz=True)
        r = self.run_script("a.fna.gz")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(len(pd.read_csv(self.d / "miss.tsv", sep="\t")), 5)

    def test_fails_when_most_loci_missing(self):
        write_fasta(self.d / "wrong.fna", ["dZ_1", "dZ_2"])
        r = self.run_script("wrong.fna")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("--max_missing_loci", r.stderr)
        # report is still written so the user can inspect it
        self.assertEqual(len(pd.read_csv(self.d / "miss.tsv", sep="\t")), 13)

    def test_threshold_is_configurable(self):
        write_fasta(self.d / "a.fna", self.present)   # 5/13 = 38% missing
        self.assertNotEqual(self.run_script("a.fna", max_frac=0.38).returncode, 0)
        self.assertEqual(self.run_script("a.fna", max_frac=0.39).returncode, 0)

    def test_unknown_dna_sample_id(self):
        write_fasta(self.d / "a.fna", self.present)
        r = self.run_script("a.fna", dna="dNOPE")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not a column", r.stderr)


if __name__ == "__main__":
    unittest.main()

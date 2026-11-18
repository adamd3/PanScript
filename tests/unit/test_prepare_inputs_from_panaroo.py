#!/usr/bin/env python3
"""Unit tests for scripts/prepare_inputs_from_panaroo.py.

The synthetic Panaroo directory follows Panaroo's own writers (panaroo/
prokka.py, find_missing.py, generate_output.py): gene_data.csv rows are
comma-joined without quoting (the description may contain commas); refound
genes are named <member>_refound_<n>; gene_presence_absence.csv cells may
carry status suffixes (_len, _pseudo, or _stop in older versions).
"""

import csv
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "prepare_inputs_from_panaroo.py"


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        lines = fh.read().splitlines()
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = ""
        else:
            seqs[name] += line
    return seqs


class PrepareInputsTest(unittest.TestCase):
    def setUp(self):
        self.d = Path(tempfile.mkdtemp(prefix="sq_panaroo_"))
        p = self.d / "panaroo"
        p.mkdir()
        long_seq = "ATG" + "GCA" * 40 + "TAA"   # > 60 nt, checks line wrapping
        rows = [
            # gff_file, scaffold, clustering_id, annotation_id, prot, dna, gene_name, description
            "isoA,ctg1,0_0_0,isoA_00001,MA,ATGGCATAA,dnaA,chromosomal replication initiator protein",
            "isoA,ctg1,0_0_1,isoA_00002,MAA,ATGGCAGCATAA,,hypothetical protein, putative",
            "isoA,ctg1,0_0_2,isoA_00003,MAAA," + long_seq + ",splX,fragment 1",
            "isoA,ctg1,0_0_3,isoA_00004,MA,ATGGCCTAA,splX,fragment 2",
            "isoA,ctg2,0_refound_7,0_refound_7,MA*,ATGTAGTAA,,location:10-18;strand:+",
            "isoB,ctg1,1_0_0,isoB_00001,MA,ATGGCGTAA,dnaA,chromosomal replication initiator protein",
            "isoB,ctg1,1_0_1,isoB_00002,MAA,ATGGCGGCATAA,,hypothetical protein",
        ]
        (p / "gene_data.csv").write_text(
            "gff_file,scaffold_name,clustering_id,annotation_id,prot_sequence,"
            "dna_sequence,gene_name,description\n" + "\n".join(rows) + "\n")
        with open(p / "gene_presence_absence.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Gene", "Non-unique Gene name", "Annotation", "isoA", "isoB"])
            w.writerow(["dnaA", "", "replication", "isoA_00001", "isoB_00001"])
            w.writerow(["group_1", "", "hypothetical", "isoA_00002_len", "isoB_00002"])
            w.writerow(["splX", "", "split gene", "isoA_00003;isoA_00004", ""])
            w.writerow(["pseudoY", "", "refound", "0_refound_7_stop", ""])
            w.writerow(["ghost", "", "no sequence", "isoA_99999_pseudo", ""])

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def run_script(self, *extra):
        return subprocess.run([sys.executable, str(SCRIPT), "--panaroo_dir", "panaroo",
                               "--outdir", "out", *extra],
                              cwd=self.d, capture_output=True, text=True)

    def test_fasta_per_isolate_with_panaroo_ids(self):
        r = self.run_script()
        self.assertEqual(r.returncode, 0, r.stderr)
        a = read_fasta(self.d / "out" / "fasta" / "isoA.fna")
        b = read_fasta(self.d / "out" / "fasta" / "isoB.fna")
        # all of isoA's genes, including the refound pseudogene, with full sequences
        self.assertEqual(set(a), {"isoA_00001", "isoA_00002", "isoA_00003", "isoA_00004",
                                  "0_refound_7"})
        self.assertEqual(a["isoA_00003"], "ATG" + "GCA" * 40 + "TAA")   # unwrapped correctly
        self.assertEqual(a["0_refound_7"], "ATGTAGTAA")
        self.assertEqual(set(b), {"isoB_00001", "isoB_00002"})
        self.assertNotIn("isoB_00001", a)

    def test_status_suffixes_matched_and_true_misses_reported(self):
        r = self.run_script()
        # isoA: 6 loci listed (suffixes stripped); 5 have sequences; ghost missing
        self.assertIn("isoA: 5 sequences written; 5/6 loci", r.stderr)
        self.assertIn("isoA_99999", r.stderr)
        self.assertIn("isoB: 2 sequences written; 2/2 loci", r.stderr)

    def test_metadata_template(self):
        self.run_script()
        with open(self.d / "out" / "metadata_template.tsv") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual([r["dna_sample_id"] for r in rows], ["isoA", "isoB"])
        self.assertTrue(rows[0]["fasta"].endswith("out/fasta/isoA.fna"))
        self.assertEqual(rows[0]["fastq1"], "")

    def test_isolate_subset_and_unknown(self):
        r = self.run_script("--isolates", "isoB")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertFalse((self.d / "out" / "fasta" / "isoA.fna").exists())
        r = self.run_script("--isolates", "nope")
        self.assertNotEqual(r.returncode, 0)


if __name__ == "__main__":
    unittest.main()

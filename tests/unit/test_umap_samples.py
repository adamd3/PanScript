#!/usr/bin/env python3
"""Unit tests for bin/umap_samples.py (PROGRESS.md P9).

Needs umap-learn and matplotlib (present in the pipeline container); skipped
otherwise. Run inside the container, e.g.:
    docker run --rm -v $PWD:/repo -w /repo <image> python -m unittest discover -s tests/unit
"""


import shutil
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "umap_samples.py"
def _umap_works():
    # a real import, not just find_spec: a broken umap/scipy install should skip
    # (with the same writable cache dirs that bin/umap_samples.py sets, so
    # this doesn't skip when run as a non-root user, as in CI)
    env = dict(os.environ)
    env.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
    env.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")
    r = subprocess.run([sys.executable, "-c", "import umap, matplotlib"],
                       capture_output=True, env=env)
    return r.returncode == 0


HAVE_UMAP = _umap_works()


@unittest.skipUnless(HAVE_UMAP, "umap-learn not installed (run inside the pipeline container)")
class UmapSamplesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = Path(tempfile.mkdtemp(prefix="sq_umap_"))
        rng = np.random.default_rng(0)
        samples = [f"S{i}" for i in range(8)]
        counts = pd.DataFrame(rng.gamma(2, 50, size=(200, 8)), columns=samples,
                              index=[f"g{i}" for i in range(200)])
        counts.iloc[:5, 0] = np.nan          # accessory genes: dropped before UMAP
        counts.index.name = "feature_id"
        counts.to_csv(cls.d / "norm.tsv", sep="\t")
        pd.DataFrame({"sample_name": samples,
                      "ST": ["1.0", "1.0", "2.0", "2.0", "3.0", "3.0", "", "x"]}
                     ).to_csv(cls.d / "meta.tsv", sep="\t", index=False)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d, ignore_errors=True)

    def run_umap(self, out, seed=42, group="ST"):
        (self.d / out).mkdir(exist_ok=True)
        return subprocess.run([sys.executable, str(SCRIPT), "--norm_counts", "norm.tsv",
                               "--metadata", "meta.tsv", "--group", group, "--seed", str(seed),
                               "--outdir", out], cwd=self.d, capture_output=True, text=True)

    def test_reproducible_with_seed(self):
        for out in ("a", "b"):
            r = self.run_umap(out)
            self.assertEqual(r.returncode, 0, r.stderr)
        a = (self.d / "a" / "umap_coordinates.tsv").read_text()
        b = (self.d / "b" / "umap_coordinates.tsv").read_text()
        self.assertEqual(a, b)
        self.assertTrue((self.d / "a" / "umap_ST.png").stat().st_size > 0)

    def test_coordinates_and_labels(self):
        self.run_umap("c")
        c = pd.read_csv(self.d / "c" / "umap_coordinates.tsv", sep="\t", dtype={"ST": str},
                        keep_default_na=False)
        self.assertEqual(list(c.columns), ["sample_name", "UMAP_1", "UMAP_2", "ST"])
        self.assertEqual(len(c), 8)
        # float-looking labels tidied; missing values shown as NA
        self.assertEqual(c.set_index("sample_name")["ST"].to_dict(),
                         {"S0": "1", "S1": "1", "S2": "2", "S3": "2", "S4": "3", "S5": "3",
                          "S6": "NA", "S7": "x"})

    def test_unknown_group_column(self):
        r = self.run_umap("d", group="nope")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not a column", r.stderr)


if __name__ == "__main__":
    unittest.main()

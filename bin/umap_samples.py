#!/usr/bin/env python3
"""UMAP of samples from size-factor-normalised expression.

Input is norm_counts.tsv (linear scale; NA = gene absent). The UMAP is run on
log2(x + 1) of the genes present in every sample (rows with any NA dropped,
i.e. the core genome), one point per sample, coloured by a metadata column.

A fixed random seed makes the embedding reproducible (umap-learn runs
single-threaded when random_state is set).

Outputs (in --outdir):
  umap_coordinates.tsv   sample_name, UMAP_1, UMAP_2, <group>
  umap_<group>.png       scatter plot
"""

import argparse
import os
import sys
import warnings

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")  # writable cache in containers
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def tidy_label(v):
    """'313.0' -> '313' (numeric metadata columns are written as floats upstream)."""
    try:
        f = float(v)
        return str(int(f)) if f.is_integer() else v
    except (TypeError, ValueError):
        return v


def n_neighbours_for(n_samples):
    # same rule as the previous R implementation, capped below n_samples
    k = 3 if n_samples < 15 else 15
    return max(2, min(k, n_samples - 1))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--norm_counts", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--group", required=True, help="metadata column used for colours")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--outdir", default=".")
    a = p.parse_args()

    counts = pd.read_csv(a.norm_counts, sep="\t", index_col=0)
    x = np.log2(counts.dropna(axis=0, how="any") + 1)
    if x.shape[1] < 3:
        sys.exit(f"ERROR: UMAP needs at least 3 samples; found {x.shape[1]}.")

    meta = pd.read_csv(a.metadata, sep="\t", dtype=str)
    if a.group not in meta.columns:
        sys.exit(f"ERROR: --group '{a.group}' is not a column of the metadata file.")
    meta = meta.drop_duplicates("sample_name").set_index("sample_name")
    groups = meta.reindex(x.columns)[a.group].fillna("NA").map(tidy_label)

    import umap  # imported late: slow (numba), and input errors should be reported first

    k = n_neighbours_for(x.shape[1])
    # compatibility_layout=True: the standard (umap-learn 0.5.x) algorithm,
    # stated explicitly so a future change of default cannot alter plots.
    # n_jobs=1 is implied by a fixed random_state (needed for reproducibility).
    reducer = umap.UMAP(n_neighbors=k, n_components=2, random_state=a.seed, n_jobs=1,
                        compatibility_layout=True)
    with warnings.catch_warnings():
        # expected with very few samples; not useful in pipeline logs
        warnings.filterwarnings("ignore", message=".*k >= N.*")
        warnings.filterwarnings("ignore", message=".*spectral initialisation failed.*")
        emb = reducer.fit_transform(x.T.to_numpy())

    coords = pd.DataFrame(emb, columns=["UMAP_1", "UMAP_2"], index=x.columns)
    coords.index.name = "sample_name"
    coords[a.group] = groups.values
    coords.to_csv(os.path.join(a.outdir, "umap_coordinates.tsv"), sep="\t")

    levels = sorted(coords[a.group].unique(), key=str)
    cmap = plt.get_cmap("tab20")
    fig, ax = plt.subplots(figsize=(9, 7))
    for i, lev in enumerate(levels):
        sub = coords[coords[a.group] == lev]
        ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=80, color=cmap(i % 20),
                   edgecolor="black", linewidth=0.8, label=str(lev))
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"Samples ({x.shape[0]} core genes, log2(norm_counts + 1), "
                 f"n_neighbors={k}, seed={a.seed})", fontsize=10)
    ax.legend(title=a.group, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, f"umap_{a.group}.png"), dpi=300)


if __name__ == "__main__":
    main()

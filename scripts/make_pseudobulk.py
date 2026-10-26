#!/usr/bin/env python3
"""Pseudobulk counts: sum UMIs over all cells of each sample (Cell Ranger filtered matrix).

Reads <cellranger_dir>/<sample>/outs/filtered_feature_bc_matrix/ (MEX format, no h5py needed) and writes
a gene x sample table in the same layout as merged_counts.csv in pipeline-bulkrnaseq-sm:
    gene_id,<sample1>,<sample2>,...
Only "Gene Expression" features are kept.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# -----------------------------
# Args
# -----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--h5", nargs="+", required=True,
                    help="<cellranger_dir>/<sample>/outs/filtered_feature_bc_matrix.h5 (the MEX folder next to it is read)")
parser.add_argument("--out", required=True)
args = parser.parse_args()

columns = {}
genes = None
ncells = {}

for h5 in args.h5:
    outs = Path(h5).parent
    sample = outs.parent.name
    mex = outs / "filtered_feature_bc_matrix"

    features = pd.read_csv(mex / "features.tsv.gz", sep="\t", header=None, usecols=[0, 1, 2],
                           names=["gene_id", "gene_symbol", "feature_type"])

    # matrix.mtx.gz: comment lines, then "n_features n_barcodes n_entries", then "feature barcode count" (1-based)
    mtx = pd.read_csv(mex / "matrix.mtx.gz", sep=" ", comment="%", header=None, names=["feature", "barcode", "count"],
                      dtype={"feature": np.int64, "barcode": np.int64, "count": np.float64})
    n_features, n_barcodes = int(mtx.iloc[0, 0]), int(mtx.iloc[0, 1])
    mtx = mtx.iloc[1:]

    if n_features != len(features):
        raise ValueError(f"{sample}: matrix has {n_features} features but features.tsv.gz has {len(features)}")

    totals = np.bincount(mtx["feature"].to_numpy() - 1, weights=mtx["count"].to_numpy(), minlength=n_features)

    gex = (features["feature_type"] == "Gene Expression").to_numpy()
    if genes is None:
        genes = features.loc[gex, "gene_id"].to_numpy()
    elif not np.array_equal(genes, features.loc[gex, "gene_id"].to_numpy()):
        raise ValueError(f"{sample}: gene list differs from previous samples (different reference?)")

    columns[sample] = totals[gex].astype(np.int64)
    ncells[sample] = n_barcodes
    print(f"[INFO] {sample}: {n_barcodes} cells, {int(totals[gex].sum()):,} UMIs")

merged = pd.DataFrame(columns, index=pd.Index(genes, name="gene_id"))
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
merged.to_csv(args.out, sep=",")

pd.Series(ncells, name="n_cells").rename_axis("Barcode").to_csv(Path(args.out).with_name("n_cells.csv"))
print(f"Pseudobulk counts: {merged.shape[0]} genes x {merged.shape[1]} samples -> {args.out}")

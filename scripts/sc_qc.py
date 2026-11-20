#!/usr/bin/env python3
"""Per-sample cell QC: QC metrics + Scrublet doublet scores for one Cell Ranger sample.

Reads <sample>/outs/filtered_feature_bc_matrix.h5 and writes one row per barcode:
    barcode, sample, n_genes, n_counts, pct_mito, doublet_score, predicted_doublet, status
status is "pass", "low_quality" (fails a QC threshold) or "doublet" (passes QC, called a doublet).
The matrix itself is not copied: sc_cluster re-reads the Cell Ranger h5 and keeps the "pass" cells.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

parser = argparse.ArgumentParser()
parser.add_argument("--h5", required=True)
parser.add_argument("--sample", required=True)
parser.add_argument("--out", required=True, help="per-cell table (.csv.gz)")
parser.add_argument("--min_genes", type=int, default=300)
parser.add_argument("--max_genes", type=int, default=10000)
parser.add_argument("--min_counts", type=int, default=500)
parser.add_argument("--max_pct_mito", type=float, default=5.0)
parser.add_argument("--doublets", default="true", help="true | false: run Scrublet")
parser.add_argument("--doublet_rate_per_1k", type=float, default=0.004,
                    help="expected doublet rate per 1,000 recovered cells (3' v4 / GEM-X ~0.004, v3.1 ~0.008)")
parser.add_argument("--seed", type=int, default=0)
args = parser.parse_args()

adata = sc.read_10x_h5(args.h5)
adata.var_names_make_unique()
n = adata.n_obs

adata.var["mt"] = adata.var_names.str.lower().str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

table = pd.DataFrame({
    "barcode": adata.obs_names,
    "sample": args.sample,
    "n_genes": adata.obs["n_genes_by_counts"].to_numpy(),
    "n_counts": adata.obs["total_counts"].to_numpy(),
    "pct_mito": adata.obs["pct_counts_mt"].to_numpy().round(3),
})

low_quality = (
    (table["n_genes"] < args.min_genes) | (table["n_genes"] > args.max_genes) |
    (table["n_counts"] < args.min_counts) | (table["pct_mito"] > args.max_pct_mito)
).to_numpy()

table["doublet_score"] = np.nan
table["predicted_doublet"] = False
if args.doublets.lower() == "true":
    # Score only cells that pass QC: empty/damaged droplets distort the simulated doublets
    ok = adata[~low_quality].copy()
    rate = min(0.25, args.doublet_rate_per_1k * ok.n_obs / 1000)
    sc.pp.scrublet(ok, expected_doublet_rate=rate, random_state=args.seed)
    table.loc[~low_quality, "doublet_score"] = ok.obs["doublet_score"].to_numpy().round(4)
    table.loc[~low_quality, "predicted_doublet"] = ok.obs["predicted_doublet"].to_numpy().astype(bool)
    print(f"[INFO] {args.sample}: Scrublet expected rate {rate:.3f}, "
          f"called {int(ok.obs['predicted_doublet'].sum())} doublets")

table["status"] = np.where(low_quality, "low_quality", np.where(table["predicted_doublet"], "doublet", "pass"))

Path(args.out).parent.mkdir(parents=True, exist_ok=True)
table.to_csv(args.out, index=False)
counts = table["status"].value_counts()
print(f"[INFO] {args.sample}: {n} cells -> pass {counts.get('pass', 0)}, "
      f"low_quality {counts.get('low_quality', 0)}, doublet {counts.get('doublet', 0)}")

#!/usr/bin/env python3
"""Summarise pseudobulk DESeq2 results across cell types and comparisons.

Reads <dge_dir>/<cell_type>/<comparison>/dge_results.csv (an empty file = comparison skipped) and writes:
    dge_summary.csv      cell_type, comparison, cells, samples, genes tested, up, down (or skip reason)
    barplot_dge_summary.png   significant genes per cell type x comparison (up / down)
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--dge_dir", required=True)
parser.add_argument("--pseudobulk_dir", required=True)
parser.add_argument("--celltypes", required=True)
parser.add_argument("--comparisons", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

types = [x for x in args.celltypes.split(",") if x]
comps = [x for x in args.comparisons.split(",") if x]
out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
ncells = pd.read_csv(Path(args.pseudobulk_dir) / "n_cells_per_sample.csv", index_col=0)

rows = []
for t in types:
    for c in comps:
        f = Path(args.dge_dir) / t / c / "dge_results.csv"
        status = (Path(args.pseudobulk_dir) / t / f"status_{c}.txt").read_text().strip()
        row = {"cell_type": t, "comparison": c,
               "cells": int(ncells.loc[t].sum()) if t in ncells.index else 0,
               "samples": len(pd.read_csv(Path(args.pseudobulk_dir) / t / "metadata.csv")),
               "status": status}
        if status == "ok" and f.stat().st_size > 0:
            d = pd.read_csv(f)
            row.update(genes_tested=int(d["padj"].notna().sum()),
                       up=int((d["direction"] == "Up").sum()), down=int((d["direction"] == "Down").sum()))
        rows.append(row)
summary = pd.DataFrame(rows)
summary.to_csv(out / "dge_summary.csv", index=False)

fig, axes = plt.subplots(1, len(comps), figsize=(3.2 * len(comps) + 3, 0.45 * len(types) + 2), squeeze=False)
for ax, c in zip(axes[0], comps):
    s = summary[summary["comparison"] == c].set_index("cell_type").reindex(types)
    up = s.get("up", pd.Series(np.nan, index=types)).to_numpy(dtype=float)
    down = s.get("down", pd.Series(np.nan, index=types)).to_numpy(dtype=float)
    y = np.arange(len(types))
    ax.barh(y, np.nan_to_num(up), color="#c0392b", height=0.7, label="up")
    ax.barh(y, -np.nan_to_num(down), color="#2c6fbb", height=0.7, label="down")
    for i in range(len(types)):
        if np.isnan(up[i]):
            ax.text(0, i, "skipped", va="center", ha="center", fontsize=7, color="#888")
        else:
            ax.text(up[i], i, f" {int(up[i])}", va="center", ha="left", fontsize=7)
            ax.text(-down[i], i, f"{int(down[i])} ", va="center", ha="right", fontsize=7)
    lim = max(1, np.nanmax(np.concatenate([up, down, [1]])) * 1.3)
    ax.set_xlim(-lim, lim); ax.axvline(0, color="#555", lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels([t.replace("_", " ") for t in types])
    ax.invert_yaxis()
    ax.set_title(c.replace("_vs_", "\nvs "), fontsize=9, loc="left")
    ax.set_xlabel("significant genes (down | up)")
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)
axes[0][0].legend(frameon=False, fontsize=8, loc="lower right")
fig.suptitle("Pseudobulk DESeq2: significant genes per cell type", x=0.01, ha="left")
fig.tight_layout(); fig.savefig(out / "barplot_dge_summary.png", dpi=150); plt.close(fig)
print(summary.to_string(index=False))

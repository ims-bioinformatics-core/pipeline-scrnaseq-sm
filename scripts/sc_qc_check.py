#!/usr/bin/env python3
"""Cell QC across samples (output of sc_qc): violins of the QC metrics and cells kept / removed per sample.

Outputs in --outdir:
    qc_summary.csv               per sample: cells in, low_quality, doublet, pass, medians of the metrics
    violin_qc_per_sample.png     genes, UMIs, % mito, doublet score per sample (thresholds dashed)
    barplot_cell_status.png      cells per sample split by pass / doublet / low_quality
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--qc", nargs="+", required=True, help="per-sample tables from sc_qc.py")
parser.add_argument("--metadata", required=True)
parser.add_argument("--group_cols", default="", help="comma-separated metadata columns used to order/colour samples")
parser.add_argument("--min_genes", type=float)
parser.add_argument("--max_genes", type=float)
parser.add_argument("--min_counts", type=float)
parser.add_argument("--max_pct_mito", type=float)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
cells = pd.concat([pd.read_csv(f) for f in args.qc], ignore_index=True)

meta = pd.read_csv(args.metadata, dtype=str)
meta.columns = [" ".join(c.split()) for c in meta.columns]
meta["Barcode"] = meta["Barcode"].str.strip()
group_cols = [c.strip() for c in args.group_cols.split(",") if c.strip() and c.strip() in meta.columns]

samples = sorted(cells["sample"].unique())
if group_cols:
    m = meta.set_index("Barcode").reindex(samples)
    m["group"] = m[group_cols].fillna("NA").agg(" | ".join, axis=1)
    samples = m.sort_values(["group"]).index.tolist()
    sample_group = m["group"].to_dict()
else:
    sample_group = {s: "all" for s in samples}
groups = list(dict.fromkeys(sample_group[s] for s in samples))
palette = plt.get_cmap("tab10")
group_colour = {g: palette(i % 10) for i, g in enumerate(groups)}

# ---- summary table ------------------------------------------------------
summary = (cells.groupby("sample")
           .agg(cells_in=("barcode", "size"),
                median_genes=("n_genes", "median"), median_counts=("n_counts", "median"),
                median_pct_mito=("pct_mito", "median")))
status = cells.groupby(["sample", "status"]).size().unstack(fill_value=0)
for s in ["pass", "doublet", "low_quality"]:
    summary[s] = status.get(s, 0)
summary["pct_pass"] = (100 * summary["pass"] / summary["cells_in"]).round(1)
summary["group"] = [sample_group[s] for s in summary.index]
summary.loc[samples].to_csv(out / "qc_summary.csv")

# ---- violins ------------------------------------------------------------
rng = np.random.default_rng(0)
plot_cells = (cells.groupby("sample", group_keys=False)
              .apply(lambda d: d.sample(min(len(d), 5000), random_state=0)))
metrics = [("n_genes", "Genes per cell", True, [args.min_genes, args.max_genes]),
           ("n_counts", "UMIs per cell", True, [args.min_counts]),
           ("pct_mito", "% mitochondrial UMIs", False, [args.max_pct_mito]),
           ("doublet_score", "Scrublet doublet score", False, [])]
fig, axes = plt.subplots(len(metrics), 1, figsize=(max(10, 0.42 * len(samples) + 3), 13), sharex=True)
for ax, (col, label, log, lines) in zip(axes, metrics):
    data = [plot_cells.loc[plot_cells["sample"] == s, col].dropna().to_numpy() for s in samples]
    data = [np.log10(d + 1) if log else d for d in data]
    pos = [i for i, d in enumerate(data) if len(d) > 1]
    parts = ax.violinplot([data[i] for i in pos], positions=pos, widths=0.8, showmedians=True, showextrema=False)
    for body, i in zip(parts["bodies"], pos):
        body.set_facecolor(group_colour[sample_group[samples[i]]]); body.set_edgecolor("none"); body.set_alpha(0.75)
    parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(1)
    for v in lines:
        if v is not None:
            ax.axhline(np.log10(v + 1) if log else v, color="#c0392b", ls="--", lw=1)
    ax.set_ylabel(("log10 " if log else "") + label)
    ax.grid(axis="y", color="#e5e5e5", lw=0.6); ax.set_axisbelow(True)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)
axes[-1].set_xticks(range(len(samples))); axes[-1].set_xticklabels(samples, rotation=90)
handles = [plt.Rectangle((0, 0), 1, 1, color=group_colour[g], alpha=0.75) for g in groups]
axes[0].legend(handles, groups, title=" | ".join(group_cols) if group_cols else None,
               loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
axes[0].set_title("Cell QC per sample (red dashed = thresholds; up to 5,000 cells per sample shown)", loc="left")
fig.tight_layout()
fig.savefig(out / "violin_qc_per_sample.png", dpi=150)
plt.close(fig)

# ---- cell status barplot -------------------------------------------------
fig, ax = plt.subplots(figsize=(max(10, 0.42 * len(samples) + 3), 5))
bottom = np.zeros(len(samples))
for s, colour in [("pass", "#3b78c2"), ("doublet", "#d08a1e"), ("low_quality", "#b0b7c0")]:
    vals = summary.loc[samples, s].to_numpy()
    ax.bar(range(len(samples)), vals, bottom=bottom, color=colour, label=s.replace("_", " "), width=0.8)
    bottom += vals
for i, smp in enumerate(samples):
    ax.text(i, bottom[i], f"{summary.loc[smp, 'pct_pass']:.0f}%", ha="center", va="bottom", fontsize=7, rotation=90)
ax.set_xticks(range(len(samples))); ax.set_xticklabels(samples, rotation=90)
for i, lab in enumerate(ax.get_xticklabels()):
    lab.set_color(group_colour[sample_group[samples[i]]])
ax.set_ylabel("Cells (Cell Ranger filtered)")
ax.set_title("Cells kept per sample (label = % passing; tick colour = group)", loc="left")
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
for side in ["top", "right"]:
    ax.spines[side].set_visible(False)
fig.tight_layout()
fig.savefig(out / "barplot_cell_status.png", dpi=150)
plt.close(fig)

tot = summary[["cells_in", "pass", "doublet", "low_quality"]].sum()
print(f"[INFO] {len(samples)} samples, {tot['cells_in']:,} cells -> pass {tot['pass']:,}, "
      f"doublet {tot['doublet']:,}, low_quality {tot['low_quality']:,}")

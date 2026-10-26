#!/usr/bin/env python3
"""Per-sample Cell Ranger QC plots from outs/metrics_summary.csv.

Outputs (in --outdir):
  cellranger_metrics.csv          one row per sample: metadata + all metrics as numbers
  barplot_cellranger_metrics.png  key metrics, one panel each, bars coloured by --group_col
  barplot_cellranger_mapping.png  read mapping composition per sample (exonic / intronic / intergenic / ...)
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# -----------------------------
# Args
# -----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--metrics", nargs="+", required=True, help="<cellranger_dir>/<sample>/outs/metrics_summary.csv files")
parser.add_argument("--metadata", required=True)
parser.add_argument("--group_col", default="", help="metadata column to colour/group samples by (empty = no grouping)")
parser.add_argument("--expected_cells", default="", help="expected cells per sample; drawn as a dashed line on Estimated Number of Cells (empty = no line)")
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)
expected_cells = float(args.expected_cells) if str(args.expected_cells).strip() else None


# -----------------------------
# Read metrics
# -----------------------------
def to_number(x):
    """'8,524' -> 8524, '93.2%' -> 93.2"""
    s = str(x).replace(",", "").replace("%", "").strip()
    try:
        return float(s)
    except ValueError:
        return np.nan


rows = []
for f in args.metrics:
    sample = Path(f).parent.parent.name          # <cellranger_dir>/<sample>/outs/metrics_summary.csv
    m = pd.read_csv(f, dtype=str).iloc[0]
    row = {"sample": sample}
    row.update({col: to_number(val) for col, val in m.items()})
    rows.append(row)

df = pd.DataFrame(rows)

# -----------------------------
# Metadata
# -----------------------------
meta = pd.read_csv(args.metadata, dtype=str, encoding="utf-8-sig")
meta.columns = meta.columns.str.strip()

if "Barcode" not in meta.columns:
    raise ValueError(f"'Barcode' column not found. Found: {list(meta.columns)}")

meta = meta.rename(columns={"Barcode": "sample"})
meta["sample"] = meta["sample"].str.strip()
df = meta.merge(df, on="sample", how="right")

# -----------------------------
# Grouping
# -----------------------------
group_col = args.group_col.strip()
if group_col:
    if group_col not in df.columns:
        raise ValueError(f"'{group_col}' not found in metadata.\nAvailable columns:\n{list(meta.columns)}")
    df["group"] = df[group_col].fillna("Unknown").replace("", "Unknown").astype(str)
else:
    df["group"] = "all"

# Label: "SITTA5 (1FH)" when a Sample name exists
if "Sample name" in df.columns:
    df["label"] = df["sample"] + df["Sample name"].fillna("").map(lambda x: f" ({x})" if x else "")
else:
    df["label"] = df["sample"]

# Order: by group, then by number of cells
df = df.sort_values(["group", "Estimated Number of Cells"], ascending=[True, False]).reset_index(drop=True)

df.drop(columns=["label", "group"]).rename(columns={"sample": "Barcode"}).to_csv(outdir / "cellranger_metrics.csv", index=False)

groups = list(dict.fromkeys(df["group"]))
palette = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#999999"]
colours = {g: palette[i % len(palette)] for i, g in enumerate(groups)}
bar_colours = df["group"].map(colours)

y = np.arange(len(df))[::-1]       # first sample on top
fig_height = max(4, len(df) * 0.28 + 1.5)

# -----------------------------
# Plot 1: key metrics
# -----------------------------
# (metric, x label, guide line or None)
panels = [
    ("Estimated Number of Cells",                    "Cells",                 expected_cells),  # expected_cells from config
    ("Mean Reads per Cell",                          "Reads",                 20000),   # 10x minimum recommendation
    ("Median Genes per Cell",                        "Genes",                 None),
    ("Median UMI Counts per Cell",                   "UMIs",                  None),
    ("Number of Reads",                              "Reads (millions)",      None),
    ("Sequencing Saturation",                        "%",                     None),
    ("Fraction Reads in Cells",                      "%",                     70),      # 10x web-summary warning threshold
    ("Valid Barcodes",                               "%",                     75),      # 10x web-summary warning threshold
    ("Reads Mapped Confidently to Genome",           "%",                     None),
    ("Reads Mapped Confidently to Transcriptome",    "%",                     30),      # 10x web-summary warning threshold
    ("Q30 Bases in RNA Read",                        "%",                     65),      # 10x web-summary warning threshold
    ("Total Genes Detected",                         "Genes",                 None),
]
panels = [p for p in panels if p[0] in df.columns]

ncol = 4
nrow = int(np.ceil(len(panels) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, fig_height * nrow * 0.85), sharey=True, squeeze=False)

for ax, (metric, xlabel, guide) in zip(axes.flat, panels):
    values = df[metric] / 1e6 if metric == "Number of Reads" else df[metric]
    ax.barh(y, values, color=bar_colours, height=0.8)
    if guide is not None:
        ax.axvline(guide, color="#555555", linestyle="--", linewidth=0.8)
    ax.set_title(metric, fontsize=10)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.tick_params(axis="x", labelsize=8)
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if xlabel == "%":
        ax.set_xlim(0, 100)

for ax in axes.flat[len(panels):]:
    ax.axis("off")

for ax in axes[:, 0]:
    ax.set_yticks(y)
    ax.set_yticklabels(df["label"], fontsize=8)

if group_col:
    handles = [plt.Rectangle((0, 0), 1, 1, color=colours[g]) for g in groups]
    fig.legend(handles, groups, title=group_col, loc="upper center", bbox_to_anchor=(0.5, 1.0),
               ncol=min(len(groups), 6), frameon=False)
fig.suptitle("Cell Ranger metrics per sample  (dashed = 10x guideline"
             + (f"; expected cells = {expected_cells:,.0f}" if expected_cells else "") + ")", y=1.02 if group_col else 1.0)
fig.tight_layout(rect=(0, 0, 1, 0.97 if group_col else 1))
fig.savefig(outdir / "barplot_cellranger_metrics.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# -----------------------------
# Plot 2: mapping composition (% of all reads)
# -----------------------------
exonic = df.get("Reads Mapped Confidently to Exonic Regions", 0)
intronic = df.get("Reads Mapped Confidently to Intronic Regions", 0)
intergenic = df.get("Reads Mapped Confidently to Intergenic Regions", 0)
genome = df.get("Reads Mapped to Genome", 0)
confident = df.get("Reads Mapped Confidently to Genome", 0)

parts = [
    ("Exonic",                 exonic,                                   "#0072B2"),
    ("Intronic",               intronic,                                 "#56B4E9"),
    ("Intergenic",             intergenic,                               "#009E73"),
    ("Mapped, not confident",  (genome - confident).clip(lower=0),       "#E69F00"),
    ("Unmapped",               (100 - genome).clip(lower=0),             "#BBBBBB"),
]

fig, ax = plt.subplots(figsize=(9, fig_height))
left = np.zeros(len(df))
for name, values, colour in parts:
    ax.barh(y, values, left=left, color=colour, label=name, height=0.8, edgecolor="white", linewidth=0.4)
    left += np.asarray(values, dtype=float)

ax.set_yticks(y)
ax.set_yticklabels(df["label"], fontsize=8)
ax.set_xlim(0, 100)
ax.set_xlabel("% of reads")
ax.set_title("Cell Ranger read mapping per sample")
for side in ("top", "right"):
    ax.spines[side].set_visible(False)

# group labels on the right, as in plot_barplot_mapping.py
if group_col:
    for g in groups:
        idx = y[df["group"] == g]
        ax.text(101, idx.mean(), g, ha="left", va="center", fontsize=9, fontweight="bold", color=colours[g])

ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=len(parts), frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(outdir / "barplot_cellranger_mapping.png", dpi=200, bbox_inches="tight")
plt.close(fig)

print(f"Cell Ranger metrics for {len(df)} samples -> {outdir}")

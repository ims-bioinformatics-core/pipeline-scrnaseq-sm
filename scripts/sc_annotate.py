#!/usr/bin/env python3
"""Annotate clusters with cell types from a marker panel, refine chosen types by subclustering,
plot cell-type composition, and write pseudobulk counts per cell type for DESeq2.

Annotation: each marker set is scored per cell (scanpy score_genes); a cluster gets the cell type with the
highest mean score (Unassigned if no score is above 0). Cells of --subcluster_types (e.g. the neurovascular
unit) are re-clustered on their own and each subcluster is labelled again, so small populations merged into
one cluster (pericytes / SMC, endothelial / fibroblasts) are separated.

Outputs in --outdir:
    cells_annotated.csv.gz          cell, sample, cluster, cell_type, UMAP coordinates (+ subcluster)
    cluster_annotation.csv          per cluster: assigned type, mean score of each type, cells
    composition.csv                 cells and proportion of each cell type per sample
    umap_celltype.png               UMAP coloured by cell type
    heatmap_cluster_scores.png      cluster x cell-type marker scores (z-scored per type)
    dotplot_celltype_markers.png    panel markers across the annotated cell types
    umap_marker_genes.png           expression of the first two markers of each type
    composition_per_sample.png      stacked cell-type proportions per sample
    composition_by_group.png        proportion of each cell type per sample, by group
    subcluster_umap.png, subcluster_dotplot.png   re-clustered --subcluster_types
In --pseudobulk_dir, per cell type in --celltypes:
    <cell_type>/counts.csv          gene_id x sample, UMIs summed over the cell type's cells
    <cell_type>/metadata.csv        metadata of the kept samples + Group (+ n_cells)
    <cell_type>/status_<comparison>.txt   "ok" or "skip: <reason>" (read by sc_dge)
    n_cells_per_sample.csv          cells per cell type and sample
"""

import argparse
import os
import re
import sys
from pathlib import Path

# Thread counts must be set before numba / BLAS start (scanpy imports numba): take --threads from argv now.
# Overrides OMP_NUM_THREADS=1 etc. inherited from a login shell.
_threads = next((sys.argv[i + 1] for i, a in enumerate(sys.argv[:-1]) if a == "--threads"), "8")
for _var in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_var] = _threads

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

parser = argparse.ArgumentParser()
parser.add_argument("--h5ad", required=True)
parser.add_argument("--markers", required=True, help="CSV: cell_type,genes (genes separated by ';')")
parser.add_argument("--metadata", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument("--pseudobulk_dir", required=True)
parser.add_argument("--plot_cols", default="")
parser.add_argument("--group_cols", default="", help="metadata columns combined into the DESeq2 Group")
parser.add_argument("--subcluster_types", default="")
parser.add_argument("--subcluster_resolution", type=float, default=0.5)
parser.add_argument("--celltypes", default="", help="cell types to write pseudobulk for")
parser.add_argument("--comparisons", default="", help="NUM_vs_REF,... (Group levels) to check replication for")
parser.add_argument("--min_cells", type=int, default=10, help="min cells of a type for a sample to enter its pseudobulk")
parser.add_argument("--threads", type=int, default=8)
parser.add_argument("--seed", type=int, default=0)
args = parser.parse_args()

sc.settings.n_jobs = args.threads
out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
pbdir = Path(args.pseudobulk_dir); pbdir.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(args.seed)
split = lambda s: [x.strip() for x in s.split(",") if x.strip()]


def label(value):
    return re.sub(r"[^A-Za-z0-9]+", "-", str(value)).strip("-")


adata = sc.read_h5ad(args.h5ad)
cluster_key = adata.uns["sc_cluster"]["cluster_key"]
use_rep = adata.uns["sc_cluster"]["use_rep"]
n_pcs = int(adata.uns["sc_cluster"]["n_pcs"])
print(f"[INFO] {adata.n_obs:,} cells, clusters from {cluster_key}")

meta = pd.read_csv(args.metadata, dtype=str)
meta.columns = [" ".join(c.split()) for c in meta.columns]
meta["Barcode"] = meta["Barcode"].str.strip()
plot_cols = [c for c in split(args.plot_cols) if c in adata.obs.columns]
group_cols = split(args.group_cols)

# -----------------------------
# Marker scores
# -----------------------------
panel = pd.read_csv(args.markers)
markers = {}
for _, row in panel.iterrows():
    genes = [g.strip() for g in str(row["genes"]).split(";") if g.strip()]
    present = [g for g in genes if g in adata.var_names]
    if len(present) < len(genes):
        print(f"[WARN] {row['cell_type']}: not in data: {sorted(set(genes) - set(present))}")
    if present:
        markers[row["cell_type"]] = present
types = list(markers)
for t, genes in markers.items():
    sc.tl.score_genes(adata, genes, score_name=f"score_{t}", random_state=args.seed)
score_cols = [f"score_{t}" for t in types]


def assign(groups):
    """Mean score per group -> best type (Unassigned if every mean is <= 0)."""
    means = adata.obs.loc[groups.index, score_cols].groupby(groups, observed=True).mean()
    means.columns = types
    best = means.idxmax(axis=1).where(means.max(axis=1) > 0, "Unassigned")
    return means, best


means, best = assign(adata.obs[cluster_key])
adata.obs["cell_type"] = adata.obs[cluster_key].map(best).astype(str)
ann = means.round(3)
ann.insert(0, "cell_type", best)
ann.insert(1, "cells", adata.obs[cluster_key].value_counts().reindex(ann.index))
ann.index.name = "cluster"
ann.to_csv(out / "cluster_annotation.csv")

# -----------------------------
# Subcluster chosen types (e.g. the neurovascular unit)
# -----------------------------
sub_types = [t for t in split(args.subcluster_types) if t in types]
adata.obs["subcluster"] = ""
if sub_types:
    sel = adata.obs["cell_type"].isin(sub_types).to_numpy()
    print(f"[INFO] subclustering {sel.sum():,} cells of {sub_types}")
    if sel.sum() >= 50:
        sub = adata[sel].copy()
        sc.pp.neighbors(sub, n_neighbors=15, n_pcs=n_pcs, use_rep=use_rep, random_state=args.seed)
        sc.tl.umap(sub, random_state=args.seed)
        sc.tl.leiden(sub, resolution=args.subcluster_resolution, key_added="subcluster", flavor="igraph",
                     n_iterations=2, directed=False, random_state=args.seed)
        sub_means = sub.obs[score_cols].groupby(sub.obs["subcluster"], observed=True).mean()
        sub_means.columns = types
        sub_best = sub_means.idxmax(axis=1).where(sub_means.max(axis=1) > 0, "Unassigned")
        sub.obs["cell_type"] = sub.obs["subcluster"].map(sub_best).astype(str)
        adata.obs.loc[sub.obs_names, "cell_type"] = sub.obs["cell_type"].to_numpy()
        adata.obs.loc[sub.obs_names, "subcluster"] = "s" + sub.obs["subcluster"].astype(str).to_numpy()
        sm = sub_means.round(3); sm.insert(0, "cell_type", sub_best); sm.index = "s" + sm.index.astype(str)
        sm.insert(1, "cells", sub.obs["subcluster"].value_counts().rename(lambda x: "s" + str(x)).reindex(sm.index))
        sm.index.name = "subcluster"
        sm.to_csv(out / "subcluster_annotation.csv")

        # plots of the subset
        cols = ["subcluster", "cell_type", "sample"] + plot_cols
        n = len(cols); ncol = 3; nrow = int(np.ceil(n / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(6.3 * ncol, 5.2 * nrow), squeeze=False)
        for ax in axes.flat[n:]:
            ax.axis("off")
        for ax, col in zip(axes.flat, cols):
            sc.pl.umap(sub, color=col, ax=ax, show=False, frameon=False, size=max(3, 30000 / sub.n_obs),
                       legend_loc="on data" if col == "subcluster" else "right margin", legend_fontsize=7,
                       title=f"{', '.join(sub_types)}: {col}" if col == "subcluster" else col)
        fig.tight_layout(); fig.savefig(out / "subcluster_umap.png", dpi=150); plt.close(fig)
        genes = {t: markers[t] for t in types if t in sub_types or t in set(sub_best)}
        dp = sc.pl.dotplot(sub, genes, groupby="subcluster", standard_scale="var", show=False, return_fig=True)
        dp.savefig(out / "subcluster_dotplot.png", dpi=150, bbox_inches="tight"); plt.close("all")
        del sub

adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
counts = adata.obs["cell_type"].value_counts()
print("[INFO] cell types:\n" + counts.to_string())

table = adata.obs[["sample", cluster_key, "subcluster", "cell_type"]].copy()
table["UMAP1"], table["UMAP2"] = adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1]
table.index.name = "cell"
table.to_csv(out / "cells_annotated.csv.gz")

# -----------------------------
# Plots
# -----------------------------
type_order = [t for t in types if t in counts.index] + [t for t in counts.index if t not in types]
pal = sc.pl.palettes.default_20
type_colour = {t: pal[i % len(pal)] for i, t in enumerate(type_order)}
adata.obs["cell_type"] = adata.obs["cell_type"].cat.reorder_categories(type_order)
adata.uns["cell_type_colors"] = [type_colour[t] for t in type_order]
umap = adata.obsm["X_umap"]
order = rng.permutation(adata.n_obs)
point = max(0.2, min(4.0, 60000 / adata.n_obs))

fig, ax = plt.subplots(figsize=(10, 8))
c = adata.obs["cell_type"].map(type_colour).astype(object).to_numpy()
ax.scatter(umap[order, 0], umap[order, 1], c=c[order], s=point, lw=0, rasterized=True)
for t in type_order:
    xy = np.median(umap[(adata.obs["cell_type"] == t).to_numpy()], axis=0)
    ax.text(xy[0], xy[1], t.replace("_", " "), fontsize=9, weight="bold", ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))
handles = [plt.Line2D([], [], marker="o", ls="", color=type_colour[t], markersize=6) for t in type_order]
ax.legend(handles, [f"{t.replace('_', ' ')} ({counts[t]:,})" for t in type_order], frameon=False,
          loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8)
ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values():
    s.set_visible(False)
ax.set_title(f"Cell types ({adata.n_obs:,} cells)", loc="left")
fig.tight_layout(); fig.savefig(out / "umap_celltype.png", dpi=150); plt.close(fig)

z = (means - means.mean()) / means.std(ddof=0).replace(0, 1)
fig, ax = plt.subplots(figsize=(0.55 * len(types) + 3, 0.3 * len(z) + 2))
im = ax.imshow(z.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
ax.set_xticks(range(len(types))); ax.set_xticklabels([t.replace("_", " ") for t in types], rotation=60, ha="right")
ax.set_yticks(range(len(z))); ax.set_yticklabels([f"{c}  →  {best[c].replace('_', ' ')}" for c in z.index], fontsize=8)
plt.colorbar(im, ax=ax, fraction=0.04, label="mean marker score (z per type)")
ax.set_title(f"Marker scores per cluster ({cluster_key}) and assigned type", loc="left")
fig.tight_layout(); fig.savefig(out / "heatmap_cluster_scores.png", dpi=150); plt.close(fig)

idx = (adata.obs.groupby("cell_type", observed=True, group_keys=False)
       .apply(lambda d: d.sample(min(len(d), 3000), random_state=args.seed)).index)
dp = sc.pl.dotplot(adata[idx], {t.replace("_", " "): g for t, g in markers.items()}, groupby="cell_type",
                   standard_scale="var", show=False, return_fig=True)
dp.savefig(out / "dotplot_celltype_markers.png", dpi=150, bbox_inches="tight"); plt.close("all")

genes = [g for t in types for g in markers[t][:2]]
ncol = 4; nrow = int(np.ceil(len(genes) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 4.0 * nrow), squeeze=False)
for ax in axes.flat[len(genes):]:
    ax.axis("off")
for ax, g in zip(axes.flat, genes):
    x = adata[:, g].X
    v = np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()
    o = order[np.argsort(v[order], kind="stable")]
    s = ax.scatter(umap[o, 0], umap[o, 1], c=v[o], s=point, lw=0, cmap="Reds",
                   vmax=max(np.percentile(v[v > 0], 99) if (v > 0).any() else 1, 0.1), rasterized=True)
    t = next(t for t in types if g in markers[t][:2])
    ax.set_title(f"{g}  ({t.replace('_', ' ')})", loc="left", fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
fig.tight_layout(); fig.savefig(out / "umap_marker_genes.png", dpi=130); plt.close(fig)

# -----------------------------
# Composition
# -----------------------------
comp = pd.crosstab(adata.obs["sample"], adata.obs["cell_type"]).reindex(columns=type_order, fill_value=0)
prop = comp.div(comp.sum(axis=1), axis=0)
m = meta.set_index("Barcode").reindex(comp.index)
group = (m[group_cols].fillna("NA").agg(" | ".join, axis=1) if group_cols
         else pd.Series("all", index=comp.index))
comp_out = comp.add_prefix("n_").join(prop.round(5).add_prefix("prop_"))
comp_out.insert(0, "group", group)
comp_out.to_csv(out / "composition.csv")

samples = group.sort_values(kind="stable").index.tolist()
fig, ax = plt.subplots(figsize=(max(10, 0.4 * len(samples) + 4), 5.5))
bottom = np.zeros(len(samples))
for t in type_order:
    v = prop.loc[samples, t].to_numpy()
    ax.bar(range(len(samples)), v, bottom=bottom, color=type_colour[t], width=0.85, label=t.replace("_", " "))
    bottom += v
ax.set_xticks(range(len(samples))); ax.set_xticklabels([f"{s}  {group[s]}" for s in samples], rotation=90, fontsize=7)
ax.set_ylabel("Proportion of cells"); ax.set_ylim(0, 1)
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8)
ax.set_title("Cell-type composition per sample", loc="left")
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(out / "composition_per_sample.png", dpi=150); plt.close(fig)

levels = sorted(group.unique())
gcol = dict(zip(levels, sc.pl.palettes.default_20))
ncol = 4; nrow = int(np.ceil(len(type_order) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(4.4 * ncol, 3.6 * nrow), squeeze=False)
for ax in axes.flat[len(type_order):]:
    ax.axis("off")
for ax, t in zip(axes.flat, type_order):
    data = [100 * prop.loc[group[group == g].index, t].to_numpy() for g in levels]
    ax.boxplot(data, widths=0.55, showfliers=False, medianprops=dict(color="black"))
    for i, (g, d) in enumerate(zip(levels, data), start=1):
        ax.scatter(i + rng.uniform(-0.15, 0.15, len(d)), d, s=16, color=gcol[g], zorder=3, edgecolor="white", lw=0.4)
    ax.set_xticks(range(1, len(levels) + 1)); ax.set_xticklabels(levels, rotation=35, ha="right", fontsize=7)
    ax.set_title(t.replace("_", " "), loc="left", fontsize=10); ax.set_ylabel("% of cells", fontsize=8)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
fig.suptitle("Cell-type proportions per sample, by group (points = samples)", x=0.01, ha="left")
fig.tight_layout(); fig.savefig(out / "composition_by_group.png", dpi=150); plt.close(fig)

# -----------------------------
# Pseudobulk per cell type
# -----------------------------
gene_ids = adata.var["gene_ids"].to_numpy() if "gene_ids" in adata.var else adata.var_names.to_numpy()
X = adata.layers["counts"]
md = meta.set_index("Barcode")
md["Group"] = (md[group_cols].apply(lambda col: col.map(label)).agg(".".join, axis=1) if group_cols else "all")
print("[INFO] Group levels: " + ", ".join(sorted(md["Group"].dropna().unique())))
comparisons = split(args.comparisons)
ncell_tab = pd.crosstab(adata.obs["cell_type"], adata.obs["sample"])
ncell_tab.to_csv(pbdir / "n_cells_per_sample.csv")

for t in split(args.celltypes):
    d = pbdir / t; d.mkdir(parents=True, exist_ok=True)
    sel = (adata.obs["cell_type"] == t).to_numpy()
    n_per = adata.obs.loc[sel, "sample"].value_counts()
    kept = [s for s in n_per.index if n_per[s] >= args.min_cells]
    cols = {}
    for s in kept:
        rows = np.where(sel & (adata.obs["sample"] == s).to_numpy())[0]
        cols[s] = np.asarray(X[rows].sum(axis=0)).ravel().round().astype(np.int64)
    pd.DataFrame(cols, index=pd.Index(gene_ids, name="gene_id")).to_csv(d / "counts.csv")
    smd = md.reindex(kept).copy()
    smd["n_cells"] = [int(n_per[s]) for s in kept]
    smd.reset_index().rename(columns={"index": "Barcode"}).to_csv(d / "metadata.csv", index=False)
    for comp_name in comparisons:
        num, ref = comp_name.split("_vs_")
        n_num, n_ref = int((smd["Group"] == num).sum()), int((smd["Group"] == ref).sum())
        status = "ok" if n_num >= 2 and n_ref >= 2 else \
            f"skip: {t} has >= {args.min_cells} cells in {n_num} {num} and {n_ref} {ref} samples (need 2 each)"
        (d / f"status_{comp_name}.txt").write_text(status + "\n")
    print(f"[INFO] pseudobulk {t}: {sel.sum():,} cells, {len(kept)} samples with >= {args.min_cells} cells")

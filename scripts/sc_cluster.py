#!/usr/bin/env python3
"""Combine QC-passed cells of all samples, integrate and cluster them, and plot the result.

Steps: keep cells with status "pass" (sc_qc) -> normalise (10k UMIs, log1p) -> highly variable genes
(per sample) -> PCA -> Harmony on --batch_key -> neighbours -> UMAP -> Leiden at several resolutions
-> marker genes per cluster (at --cluster_resolution).

Outputs in --outdir:
    adata_clustered.h5ad          X = log-normalised, layers["counts"] = raw UMIs, obs = metadata + QC + clusters
    cluster_markers.csv           top marker genes per cluster (Wilcoxon, up to 1,000 cells per cluster)
    pca_pseudobulk_samples.csv    sample-level PCA of pseudobulk counts
    pca_variance.png              variance explained per PC (dashed = PCs used)
    pca_cells.png                 PC1 vs PC2 of cells before Harmony, by sample and metadata
    pca_pseudobulk_samples.png    sample-level PCA (like the bulk pipeline), by metadata
    umap_clusters.png             UMAP coloured by cluster at each resolution
    umap_metadata.png             UMAP coloured by sample and each --plot_cols column
    umap_qc.png                   UMAP coloured by genes, UMIs, % mito, doublet score
    umap_split_groups.png         one panel per --split_cols group (cells of that group in colour)
    barplot_cluster_samples.png   cluster sizes and the share of each sample per cluster (integration check)
    dotplot_cluster_markers.png   top 3 marker genes per cluster
"""

import argparse
import os
import sys

# Thread counts must be set before numba / BLAS start (scanpy imports numba): take --threads from argv now.
# Overrides OMP_NUM_THREADS=1 etc. inherited from a login shell.
_threads = next((sys.argv[i + 1] for i, a in enumerate(sys.argv[:-1]) if a == "--threads"), "8")
for _var in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_var] = _threads
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

parser = argparse.ArgumentParser()
parser.add_argument("--h5", nargs="+", required=True)
parser.add_argument("--qc", nargs="+", required=True)
parser.add_argument("--metadata", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument("--n_hvg", type=int, default=3000)
parser.add_argument("--n_pcs", type=int, default=30)
parser.add_argument("--batch_key", default="sample", help="obs column for Harmony ('sample' or a metadata column; 'none' = no integration)")
parser.add_argument("--n_neighbors", type=int, default=15)
parser.add_argument("--resolutions", default="0.2,0.5,1.0")
parser.add_argument("--cluster_resolution", type=float, default=0.5)
parser.add_argument("--plot_cols", default="")
parser.add_argument("--split_cols", default="", help="metadata columns combined into groups for umap_split_groups.png")
parser.add_argument("--subsample", default="", help="optional: keep at most N cells per sample (for testing)")
parser.add_argument("--threads", type=int, default=8)
parser.add_argument("--seed", type=int, default=0)
args = parser.parse_args()

sc.settings.n_jobs = args.threads
out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(args.seed)
resolutions = [float(r) for r in args.resolutions.split(",") if r.strip()]
if args.cluster_resolution not in resolutions:
    resolutions.append(args.cluster_resolution)
cluster_key = f"leiden_{args.cluster_resolution:g}"


def clean_cols(text, columns):
    cols = [c.strip() for c in text.split(",") if c.strip()]
    missing = [c for c in cols if c not in columns]
    if missing:
        raise ValueError(f"Metadata column(s) not found: {missing}")
    return cols


# -----------------------------
# Load QC-passed cells
# -----------------------------
meta = pd.read_csv(args.metadata, dtype=str)
meta.columns = [" ".join(c.split()) for c in meta.columns]
meta["Barcode"] = meta["Barcode"].str.strip()
meta = meta.set_index("Barcode")
plot_cols = clean_cols(args.plot_cols, meta.columns)
split_cols = clean_cols(args.split_cols, meta.columns)

qc = {Path(f).name.replace("_cells.csv.gz", ""): f for f in args.qc}
parts = []
for h5 in args.h5:
    sample = Path(h5).parent.parent.name
    a = sc.read_10x_h5(h5)
    a.var_names_make_unique()
    cells = pd.read_csv(qc[sample]).set_index("barcode")
    keep = cells.index[cells["status"] == "pass"]
    if args.subsample:
        n = int(args.subsample)
        if len(keep) > n:
            keep = pd.Index(rng.choice(keep, n, replace=False))
    a = a[keep].copy()
    a.obs = cells.loc[keep, ["n_genes", "n_counts", "pct_mito", "doublet_score"]].copy()
    a.obs["sample"] = sample
    for col in meta.columns:
        a.obs[col] = meta.loc[sample, col] if sample in meta.index else "NA"
    a.obs_names = [f"{sample}_{b}" for b in a.obs_names]
    a.X = a.X.astype(np.float32)
    parts.append(a)
    print(f"[INFO] {sample}: {a.n_obs} cells")

adata = ad.concat(parts, join="inner", merge="same")
del parts
for col in adata.obs.columns:
    if adata.obs[col].dtype == object:
        adata.obs[col] = adata.obs[col].astype("category")
print(f"[INFO] combined: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

# -----------------------------
# Normalise, HVG, PCA, Harmony, UMAP, Leiden
# -----------------------------
sc.pp.filter_genes(adata, min_cells=10)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

batch_key = None if args.batch_key.lower() == "none" else args.batch_key
sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, flavor="seurat",
                            batch_key="sample" if adata.obs["sample"].nunique() > 1 else None)
hv = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(hv, max_value=10)
sc.tl.pca(hv, n_comps=max(50, args.n_pcs), svd_solver="arpack", random_state=args.seed)
adata.obsm["X_pca"] = hv.obsm["X_pca"].astype(np.float32)
adata.uns["pca"] = hv.uns["pca"]
del hv

rep = "X_pca"
if batch_key:
    sc.external.pp.harmony_integrate(adata, key=batch_key, basis="X_pca", adjusted_basis="X_pca_harmony",
                                     max_iter_harmony=20, random_state=args.seed)
    rep = "X_pca_harmony"
    print(f"[INFO] Harmony on '{batch_key}'")

sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs, use_rep=rep, random_state=args.seed)
sc.tl.umap(adata, random_state=args.seed)
for r in resolutions:
    key = f"leiden_{r:g}"
    sc.tl.leiden(adata, resolution=r, key_added=key, flavor="igraph", n_iterations=2,
                 directed=False, random_state=args.seed)
    print(f"[INFO] {key}: {adata.obs[key].nunique()} clusters")

# -----------------------------
# Plot helpers
# -----------------------------
order = rng.permutation(adata.n_obs)
umap = adata.obsm["X_umap"]
point = max(0.2, min(4.0, 60000 / adata.n_obs))


def palette(n):
    base = sc.pl.palettes.default_20 if n <= 20 else sc.pl.palettes.default_102
    return [base[i % len(base)] for i in range(n)]


def style(ax, title):
    ax.set_title(title, loc="left", fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
    for side in ax.spines.values():
        side.set_visible(False)


def umap_cat(ax, col, title=None, labels_on_data=False, legend=True):
    cats = adata.obs[col].astype("category")
    levels = list(cats.cat.categories)
    colours = dict(zip(levels, palette(len(levels))))
    c = cats.map(colours).astype(object).to_numpy()
    ax.scatter(umap[order, 0], umap[order, 1], c=c[order], s=point, lw=0, rasterized=True)
    if labels_on_data:
        for lev in levels:
            xy = np.median(umap[(cats == lev).to_numpy()], axis=0)
            ax.text(xy[0], xy[1], lev, fontsize=8, weight="bold", ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.6))
    elif legend:
        handles = [plt.Line2D([], [], marker="o", ls="", color=colours[l], markersize=5) for l in levels]
        ax.legend(handles, levels, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False,
                  fontsize=7, ncol=1 if len(levels) <= 18 else 2, markerscale=1)
    style(ax, title or col)


def umap_num(ax, values, title, cmap="viridis"):
    v = np.asarray(values, dtype=float)
    lo, hi = np.nanpercentile(v, [1, 99])
    o = order[np.argsort(v[order], kind="stable")]  # high values on top
    s = ax.scatter(umap[o, 0], umap[o, 1], c=v[o], s=point, lw=0, cmap=cmap, vmin=lo, vmax=hi, rasterized=True)
    plt.colorbar(s, ax=ax, fraction=0.04, pad=0.01)
    style(ax, title)


def grid(n, ncol=3, w=6.2, h=5.0):
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(w * ncol, h * nrow), squeeze=False)
    for ax in axes.flat[n:]:
        ax.axis("off")
    return fig, axes.flat


# -----------------------------
# PCA plots
# -----------------------------
var = adata.uns["pca"]["variance_ratio"]
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(np.arange(1, len(var) + 1), 100 * var, "o-", ms=3, color="#3b78c2")
ax.axvline(args.n_pcs, color="#c0392b", ls="--", lw=1)
ax.set_xlabel("Principal component"); ax.set_ylabel("% variance explained")
ax.set_title(f"PCA of {args.n_hvg} variable genes (dashed = {args.n_pcs} PCs used)", loc="left")
for side in ["top", "right"]:
    ax.spines[side].set_visible(False)
fig.tight_layout(); fig.savefig(out / "pca_variance.png", dpi=150); plt.close(fig)

pca_cols = ["sample"] + plot_cols
fig, axes = grid(len(pca_cols))
pc = adata.obsm["X_pca"]
for ax, col in zip(axes, pca_cols):
    cats = adata.obs[col].astype("category")
    colours = dict(zip(cats.cat.categories, palette(len(cats.cat.categories))))
    ax.scatter(pc[order, 0], pc[order, 1], c=cats.map(colours).astype(object).to_numpy()[order],
               s=point, lw=0, rasterized=True)
    handles = [plt.Line2D([], [], marker="o", ls="", color=colours[l], markersize=5) for l in cats.cat.categories]
    ax.legend(handles, cats.cat.categories, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, fontsize=7,
              ncol=1 if len(colours) <= 18 else 2)
    ax.set_xlabel(f"PC1 ({100 * var[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({100 * var[1]:.1f}%)")
    ax.set_title(f"{col} (cells, before Harmony)", loc="left", fontsize=11)
fig.tight_layout(); fig.savefig(out / "pca_cells.png", dpi=130); plt.close(fig)

# Sample-level PCA of pseudobulk counts (like the bulk pipeline)
samples = list(adata.obs["sample"].cat.categories)
onehot = sp.csr_matrix((np.ones(adata.n_obs), (adata.obs["sample"].cat.codes.to_numpy(), np.arange(adata.n_obs))),
                       shape=(len(samples), adata.n_obs))
pb = np.asarray((onehot @ adata.layers["counts"]).todense())
cpm = np.log2(pb / pb.sum(axis=1, keepdims=True) * 1e6 + 1)
top = np.argsort(cpm.var(axis=0))[::-1][:2000]
x = cpm[:, top] - cpm[:, top].mean(axis=0)
u, s, _ = np.linalg.svd(x, full_matrices=False)
pcs = u * s
pvar = s ** 2 / np.sum(s ** 2)
pb_pca = pd.DataFrame(pcs[:, :min(10, pcs.shape[1])], index=samples,
                      columns=[f"PC{i + 1}" for i in range(min(10, pcs.shape[1]))])
pb_pca = pb_pca.join(meta, how="left")
pb_pca.index.name = "Barcode"
pb_pca.to_csv(out / "pca_pseudobulk_samples.csv")
cols = plot_cols or ["sample"]
fig, axes = grid(len(cols), ncol=min(3, len(cols)), w=6.0, h=5.0)
for ax, col in zip(axes, cols):
    vals = pb_pca[col].fillna("NA") if col in pb_pca.columns else pd.Series(samples, index=samples)
    levels = sorted(vals.unique(), key=str)
    colours = dict(zip(levels, palette(len(levels))))
    for lev in levels:
        sel = (vals == lev).to_numpy()
        ax.scatter(pb_pca["PC1"][sel], pb_pca["PC2"][sel], s=60, color=colours[lev], label=lev, edgecolor="white", lw=0.6)
    for smp, row in pb_pca.iterrows():
        ax.annotate(smp, (row["PC1"], row["PC2"]), fontsize=6, xytext=(3, 3), textcoords="offset points", color="#555")
    ax.set_xlabel(f"PC1 ({100 * pvar[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({100 * pvar[1]:.1f}%)")
    ax.legend(frameon=False, fontsize=8, loc="best"); ax.set_title(f"Pseudobulk samples: {col}", loc="left", fontsize=11)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)
fig.tight_layout(); fig.savefig(out / "pca_pseudobulk_samples.png", dpi=150); plt.close(fig)

# -----------------------------
# UMAP plots
# -----------------------------
fig, axes = grid(len(resolutions), ncol=min(3, len(resolutions)), w=6.5, h=6.0)
for ax, r in zip(axes, resolutions):
    key = f"leiden_{r:g}"
    umap_cat(ax, key, f"Leiden resolution {r:g} ({adata.obs[key].nunique()} clusters)", labels_on_data=True)
fig.tight_layout(); fig.savefig(out / "umap_clusters.png", dpi=150); plt.close(fig)

meta_cols = ["sample"] + plot_cols
fig, axes = grid(len(meta_cols))
for ax, col in zip(axes, meta_cols):
    umap_cat(ax, col)
fig.tight_layout(); fig.savefig(out / "umap_metadata.png", dpi=150); plt.close(fig)

fig, axes = grid(4, ncol=2, w=6.5, h=5.2)
for ax, (col, title) in zip(axes, [("n_genes", "Genes per cell"), ("n_counts", "UMIs per cell"),
                                   ("pct_mito", "% mitochondrial UMIs"), ("doublet_score", "Scrublet doublet score")]):
    umap_num(ax, adata.obs[col], title)
fig.tight_layout(); fig.savefig(out / "umap_qc.png", dpi=150); plt.close(fig)

if split_cols:
    grp = adata.obs[split_cols].astype(str).agg(" | ".join, axis=1)
    levels = sorted(grp.unique())
    fig, axes = grid(len(levels), ncol=min(4, len(levels)), w=5.0, h=4.8)
    for ax, lev, colour in zip(axes, levels, palette(len(levels))):
        sel = (grp == lev).to_numpy()
        ax.scatter(umap[~sel, 0], umap[~sel, 1], c="#e3e6ea", s=point, lw=0, rasterized=True)
        o = order[sel[order]]
        ax.scatter(umap[o, 0], umap[o, 1], c=colour, s=point, lw=0, rasterized=True)
        style(ax, f"{lev}  (n = {sel.sum():,})")
    fig.tight_layout(); fig.savefig(out / "umap_split_groups.png", dpi=150); plt.close(fig)

# Cluster sizes + sample share per cluster
tab = pd.crosstab(adata.obs[cluster_key], adata.obs["sample"])
share = tab.div(tab.sum(axis=1), axis=0)
fig, (a1, a2) = plt.subplots(2, 1, figsize=(max(9, 0.35 * len(tab) + 4), 8), sharex=True,
                             gridspec_kw={"height_ratios": [1, 2]})
a1.bar(range(len(tab)), tab.sum(axis=1), color="#3b78c2")
a1.set_ylabel("Cells"); a1.set_yscale("log")
a1.set_title(f"Clusters ({cluster_key}): size and share of each sample", loc="left")
bottom = np.zeros(len(tab))
for smp, colour in zip(share.columns, palette(len(share.columns))):
    a2.bar(range(len(tab)), share[smp], bottom=bottom, color=colour, width=0.85, label=smp)
    bottom += share[smp].to_numpy()
a2.set_ylabel("Share of cells"); a2.set_xticks(range(len(tab))); a2.set_xticklabels(tab.index)
a2.set_xlabel("Cluster")
a2.legend(frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(1.0, 1.0), ncol=2)
for ax in (a1, a2):
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)
fig.tight_layout(); fig.savefig(out / "barplot_cluster_samples.png", dpi=150); plt.close(fig)

# -----------------------------
# Cluster markers (subsampled for speed)
# -----------------------------
idx = (adata.obs.groupby(cluster_key, observed=True, group_keys=False)
       .apply(lambda d: d.sample(min(len(d), 1000), random_state=args.seed)).index)
sub = adata[idx].copy()
sc.tl.rank_genes_groups(sub, groupby=cluster_key, method="wilcoxon", pts=True, use_raw=False)
markers = sc.get.rank_genes_groups_df(sub, group=None)
markers = markers[(markers["logfoldchanges"] > 0)].groupby("group", observed=True).head(50)
markers = markers.rename(columns={"group": "cluster", "names": "gene"})
markers.to_csv(out / "cluster_markers.csv", index=False)
top3 = {c: markers.loc[markers["cluster"] == c, "gene"].head(3).tolist() for c in tab.index}
top3 = {c: g for c, g in top3.items() if g}
dp = sc.pl.dotplot(sub, top3, groupby=cluster_key, standard_scale="var", show=False, return_fig=True)
dp.savefig(out / "dotplot_cluster_markers.png", dpi=150, bbox_inches="tight")
plt.close("all")
del sub

# -----------------------------
# Save
# -----------------------------
adata.obs["sample"] = adata.obs["sample"].astype("category")
adata.uns.pop("neighbors", None)
for k in ["distances", "connectivities"]:
    adata.obsp.pop(k, None)
adata.uns["sc_cluster"] = {"cluster_key": cluster_key, "batch_key": batch_key or "none", "n_pcs": args.n_pcs,
                           "use_rep": rep}
adata.write_h5ad(out / "adata_clustered.h5ad", compression="lzf")
print(f"[INFO] saved {out / 'adata_clustered.h5ad'}")

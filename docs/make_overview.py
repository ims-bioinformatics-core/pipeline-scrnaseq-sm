"""Draw the pipeline overview. Usage: python docs/make_overview.py docs/pipeline_overview.pdf  (or .svg / .png)

Same style as pipeline-bulkrnaseq-sm/docs/make_overview.py.
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

out = sys.argv[1]
plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})

INK = "#1f2933"; MUTED = "#5f6b7a"
C = {  # fill, edge
    "input":    ("#eef1f5", "#8a96a6"),
    "process":  ("#e3eefb", "#3b78c2"),
    "qc":       ("#fdf1df", "#d08a1e"),
    "check":    ("#f1ebfb", "#7b55c7"),
    "analysis": ("#e4f4ea", "#2f8f57"),
}
GATE = C["qc"][1]

fig = plt.figure(figsize=(11.69, 8.27))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 70.7); ax.axis("off")


def box(x, y, w, h, kind, title, tool, lines, title_fs=10.5, small=False):
    f, e = C[kind]
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0,rounding_size=1.0",
                                fc=f, ec=e, lw=1.4))
    top = y + h / 2
    d = 0.85 if small else 1.0
    ax.text(x, top - 1.9 * d, title, ha="center", va="center", fontsize=title_fs, weight="bold", color=INK)
    ax.text(x, top - 4.0 * d, tool, ha="center", va="center", fontsize=8.3, color=e, weight="bold")
    ax.text(x, top - 5.5 * d, "\n".join(lines), ha="center", va="top", fontsize=7.4, color=MUTED, linespacing=1.4)


def arrow(p, q, dashed=False, color=INK, rad=0.0):
    ax.annotate("", xy=q, xytext=p, arrowprops=dict(
        arrowstyle="-|>", lw=1.3, color=color, linestyle=(0, (3, 2)) if dashed else "-",
        shrinkA=0, shrinkB=0, mutation_scale=11, connectionstyle=f"arc3,rad={rad}"))


def elbow(points, color=INK):
    """Right-angled connector through points; arrowhead on the last segment."""
    xs_, ys_ = zip(*points[:-1])
    ax.plot(xs_, ys_, color=color, lw=1.3, solid_joinstyle="miter", solid_capstyle="butt")
    arrow(points[-2], points[-1], color=color)


# --- title ---------------------------------------------------------------
ax.text(4, 67.4, "pipeline-scrnaseq-sm", fontsize=20, weight="bold", color=INK, va="center")
ax.text(4, 64.2, "10x Genomics single-cell RNA-seq on SLURM with Snakemake: FASTQ to cell types, "
        "pseudobulk differential expression and pathways", fontsize=10.5, color=MUTED, va="center")

# --- stage bands ---------------------------------------------------------
bands = [(1.75, 14.25, "Inputs"), (15.75, 28.25, "1  Read QC + mapping"), (29.75, 42.25, "2  Cell QC"),
         (43.75, 70.25, "3  Clustering + cell types"), (71.75, 98.25, "4  Pseudobulk DE")]
for x0, x1, label in bands:
    ax.add_patch(FancyBboxPatch((x0, 57.4), x1 - x0, 3.4, boxstyle="round,pad=0,rounding_size=0.8",
                                fc="#f4f6f8", ec="none"))
    ax.text((x0 + x1) / 2, 59.1, label, ha="center", va="center", fontsize=8.6, weight="bold", color=INK)

# --- main path -----------------------------------------------------------
Y, H, W = 46.5, 12.5, 12.5
xs = [8, 22, 36, 50, 64, 78, 92]
box(xs[0], Y, W, H, "input", "Inputs", "you provide",
    ["FASTQ R1+R2 per lane", "inputs/filepaths.csv", "inputs/metadata.csv", "config.yml"])
box(xs[1], Y, W, H, "process", "cellranger_count", "Cell Ranger", ["all lanes of a sample", "UMI counts per cell", "transcriptome"],
    title_fs=9.6)
box(xs[2], Y, W, H, "process", "sc_qc", "scanpy · Scrublet", ["genes, UMIs, % mito", "doublet detection", "keep / drop per cell"])
box(xs[3], Y, W, H, "analysis", "sc_cluster", "scanpy · Harmony", ["HVGs, PCA, Harmony", "UMAP, Leiden", "sc_resolutions"])
box(xs[4], Y, W, H, "analysis", "sc_annotate", "marker scores", ["cell type per cluster", "sc_markers CSV", "pseudobulk per type"])
box(xs[5], Y, W, H, "analysis", "sc_dge", "DESeq2", ["pseudobulk per type", "sc_dge_comparisons", "results, volcano"])
box(xs[6], Y, W, H, "analysis", "sc_pathway", "clusterProfiler", ["GSEA per comparison", "GO · KEGG · Reactome", "dotplot"])
for a, b in zip(xs[:-1], xs[1:]):
    arrow((a + W / 2, Y), (b - W / 2, Y))
MB = Y - H / 2  # bottom of main row

# --- lower rows ----------------------------------------------------------
R1, R2, R3, QH = 33.5, 23.0, 12.5, 8.6
QW = 10.6
top = lambda y: y + QH / 2
bot = lambda y: y - QH / 2

# read QC: Inputs -> fastqc -> multiqc -> (gate) cellranger; multiqc -> qc_check
fq_x, mq_x = 8, 20.6
box(fq_x, R1, QW, QH, "qc", "fastqc", "FastQC", ["each FASTQ pair"], 9.5, small=True)
box(mq_x, R1, QW, QH, "qc", "multiqc", "MultiQC", ["combined report"], 9.5, small=True)
box(mq_x, R2, QW, QH, "check", "qc_check", "R", ["by run + lane,", "condition, sex"], 9.5, small=True)
arrow((fq_x, MB), (fq_x, top(R1)))
arrow((fq_x + QW / 2, R1), (mq_x - QW / 2, R1))
arrow((mq_x, top(R1)), (mq_x, MB), dashed=True, color=GATE)
arrow((mq_x, bot(R1)), (mq_x, top(R2)))

# Cell Ranger outputs: cellranger_check, merge_counts -> sex_check
cr_x = 33.6
drop = 27.2
box(cr_x, R1, QW, QH, "check", "cellranger_check", "Python", ["cells, reads,", "mapping per sample"], 8.6, small=True)
box(cr_x, R2, QW, QH, "process", "merge_counts", "Python", ["pseudobulk", "genes × samples"], 9.5, small=True)
box(cr_x, R3, QW, QH, "check", "sex_check", "R", ["Xist vs Y genes", "vs metadata sex"], 9.5, small=True)
ax.plot([drop, drop], [MB, R2], color=INK, lw=1.3, solid_capstyle="butt")
arrow((drop, R1), (cr_x - QW / 2, R1))
arrow((drop, R2), (cr_x - QW / 2, R2))
arrow((cr_x, bot(R2)), (cr_x, top(R3)))

# cell QC summary gates clustering
qc_x = 45.6
box(qc_x, R1, QW, QH, "qc", "sc_qc_check", "Python", ["cells kept per", "sample, violins"], 9.5, small=True)
arrow((xs[2] + 3, MB), (qc_x - 2.5, top(R1)))
arrow((qc_x + 2.5, top(R1)), (xs[3] - 1.5, MB), dashed=True, color=GATE)

# cluster / annotate plots (outputs of those two rules)
px, pw, ph = 65.5, 18.5, 17
py = top(R1) - ph / 2
box(px, py, pw, ph, "check", "", "", [], 9.5)
ax.text(px, py + ph / 2 - 1.9, "figures", ha="center", va="center", fontsize=9.5, weight="bold", color=INK)
ax.text(px, py + ph / 2 - 4.0, "sc_cluster + sc_annotate", ha="center", va="center", fontsize=8.3,
        weight="bold", color=C["check"][1])
ax.text(px, py + ph / 2 - 5.5, "\n".join([
    "UMAP: clusters, QC, metadata",
    "PCA: cells, pseudobulk samples",
    "cluster marker genes",
    "UMAP of cell types",
    "marker-gene dotplot",
    "cell-type composition",
]), ha="center", va="top", fontsize=7.4, color=MUTED, linespacing=1.45)
arrow((xs[3] + 3, MB), (px - pw / 2 + 1.5, py + ph / 2))
arrow((xs[4], MB), (xs[4], py + ph / 2))

# DE summary
box(xs[5] + 7, R1, 13, QH, "check", "sc_dge_summary", "Python", ["genes up / down", "per type × comparison"], 9.5, small=True)
arrow((xs[5], MB), (xs[5] + 2, top(R1)))

# --- scope brackets ------------------------------------------------------
def bracket(x0, x1, y, label):
    ax.plot([x0, x0, x1, x1], [y - 0.8, y, y, y - 0.8], color=MUTED, lw=1)
    ax.text((x0 + x1) / 2, y + 0.9, label, ha="center", va="bottom", fontsize=7.6, color=MUTED, style="italic")

by = Y + H / 2 + 1.2
bracket(xs[1] - W / 2, xs[2] + W / 2, by, "once per sample")
bracket(xs[3] - W / 2, xs[4] + W / 2, by, "all samples together (integrated)")
bracket(xs[5] - W / 2, xs[6] + W / 2, by, "per cell type × comparison")

# --- legend + footer -----------------------------------------------------
ax.add_patch(FancyBboxPatch((4, 6.4), 92, 0.08, boxstyle="square,pad=0", fc="#d5dbe2", ec="none"))
ly = 4.0
x = 4
for k, label in [("input", "input"), ("process", "processing"), ("qc", "QC report"), ("check", "check plots"), ("analysis", "analysis")]:
    f, e = C[k]
    ax.add_patch(FancyBboxPatch((x, ly - 1), 2.6, 2, boxstyle="round,pad=0,rounding_size=0.4", fc=f, ec=e, lw=1.2))
    ax.text(x + 3.4, ly, label, va="center", fontsize=8.5, color=INK)
    x += 11.8
arrow((x, ly), (x + 5, ly)); ax.text(x + 6, ly, "data flow", va="center", fontsize=8.5, color=INK)
x += 13.5
arrow((x, ly), (x + 5, ly), dashed=True, color=GATE)
ax.text(x + 6, ly, "gate: next step waits for this", va="center", fontsize=8.5, color=INK)

ax.text(4, 1.3, "Run:  sbatch submit_slurm.sh   (dry run: bash run.sh -n)    Results:  <run_dir>/<dir_flag>/    "
        "Test data: 10x E18 mouse brain (neuron_1k_v3) split into 4 samples",
        fontsize=8, color=MUTED, va="center", family="DejaVu Sans Mono")

fig.savefig(out)
print("wrote", out)

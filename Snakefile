import os
import re


# -----------------------------
# Run directory
# -----------------------------
RUN_DIR = config.get("run_dir", ".")  # default is current directory
SCRIPTS_DIR = os.path.join(config.get("pipeline_dir"), "scripts")
DIR_FLAG = config.get("dir_flag", "")


# -----------------------------
# Exclude samples (optional)
# -----------------------------
exclude_raw = config.get("exclude", "")
exclude_samples = set(
    x.strip() for x in exclude_raw.split(",") if x.strip()
)

print("[INFO] Excluding samples:", exclude_samples)


# -----------------------------
# Read filepaths and metadata
# -----------------------------
filepaths = config.get("filepaths")
if filepaths is None:
    raise ValueError(
        "Config key 'filepaths' not found. Make sure you are passing the correct config.yml with --configfile."
    )
filepaths = os.path.join(RUN_DIR, filepaths)

# filepaths.csv: sample,fastq_1,fastq_2 - one row per FASTQ pair.
# A sample sequenced on several lanes / runs has several rows; all its pairs go into one `cellranger count`.
samples = {}
with open(filepaths) as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split(",")
        if len(parts) < 3:
            continue

        sample = parts[0]

        # APPLY EXCLUSION HERE
        if sample in exclude_samples:
            print(f"[INFO] Skipping excluded sample: {sample}")
            continue

        # Cell Ranger only accepts letters, digits, '-' and '_' in --id / --sample
        if not re.match(r"^[A-Za-z0-9_-]+$", sample):
            raise ValueError(f"Sample name '{sample}' is not valid for Cell Ranger (allowed: letters, digits, '-', '_').")

        fq1 = os.path.join(RUN_DIR, parts[1])
        fq2 = os.path.join(RUN_DIR, parts[2])
        if any(p["r1"] == fq1 for p in samples.get(sample, [])):
            raise ValueError(f"Duplicate row for sample '{sample}' in {filepaths}: {parts[1]}")
        samples.setdefault(sample, []).append({"r1": fq1, "r2": fq2})

# One FastQC "unit" per FASTQ pair: <sample> if the sample has one pair, else <sample>_L001, <sample>_L002, ...
units = {}
for sample, pairs in samples.items():
    for i, pair in enumerate(pairs, start=1):
        unit = sample if len(pairs) == 1 else f"{sample}_L{i:03d}"
        units[unit] = pair

print("[INFO] Samples:", ", ".join(f"{s} ({len(p)})" if len(p) > 1 else s for s, p in samples.items()))
print(f"[INFO] {len(samples)} samples, {len(units)} FASTQ pairs")

metadata = os.path.join(RUN_DIR, config["metadata"])
# -----------------------------
# Directories for outputs under RUN_DIR
# -----------------------------

DIR_FLAG = config["dir_flag"]

fastqc_dir        = os.path.join(RUN_DIR, DIR_FLAG, "fastq-raw")
multiqc_dir       = os.path.join(RUN_DIR, DIR_FLAG, "multiqc-raw")
qc_check_dir      = os.path.join(RUN_DIR, DIR_FLAG, "qc_check")
cellranger_dir    = os.path.join(RUN_DIR, DIR_FLAG, "cellranger")
cellranger_check_dir = os.path.join(RUN_DIR, DIR_FLAG, "cellranger_check")
merged_counts_dir = os.path.join(RUN_DIR, DIR_FLAG, "merged-counts")
sample_check_dir  = os.path.join(RUN_DIR, DIR_FLAG, "sample_check")

# -----------------------------
# Parameters
# -----------------------------
fastqc_threads          = int(config["fastqc_threads"])
multiqc_mode            = config["multiqc_mode"]
cellranger_bin          = config["cellranger_bin"]
cellranger_ref          = config["cellranger_ref"]
cellranger_threads      = int(config["cellranger_threads"])
cellranger_mem          = int(config["cellranger_mem"])        # GB
cellranger_chemistry    = config["cellranger_chemistry"]
cellranger_create_bam   = str(config["cellranger_create_bam"]).lower()
cellranger_introns      = str(config["cellranger_include_introns"]).lower()
cellranger_expect_cells = config.get("cellranger_expect_cells", "") or ""
cellranger_check_cat    = config.get("cellranger_check_cat", "") or ""   # metadata column to group plots by
expected_cells          = config.get("expected_cells", "") or ""         # plot only: dashed line on Estimated Number of Cells
annotation              = config["annotation"]                           # gene annotation TSV (ensembl_gene_id, gene_symbol)
qc_check_cat            = [x.strip() for x in (config.get("qc_check_cat", "") or "").split(",") if x.strip()]   # metadata columns for the FastQC plots

# -----------------------------
# Downstream single-cell analysis (sc_* rules)
# -----------------------------
# Results go to <run_dir>/<dir_flag>/<sc_analysis>/ so several parameter sets can sit side by side
# without touching the Cell Ranger outputs.
def _csv(key, default=""):
    v = config.get(key, default)
    v = ",".join(map(str, v)) if isinstance(v, (list, tuple)) else ("" if v is None else str(v))
    return [x.strip() for x in v.split(",") if x.strip()]

sc_downstream     = str(config.get("sc_downstream", True)).lower() == "true"
sc_dir            = os.path.join(RUN_DIR, DIR_FLAG, config.get("sc_analysis", "downstream_analysis"))
sc_qc_dir         = os.path.join(sc_dir, "sc_qc")
sc_qc_check_dir   = os.path.join(sc_dir, "sc_qc_check")
sc_cluster_dir    = os.path.join(sc_dir, "sc_cluster")
sc_annotate_dir   = os.path.join(sc_dir, "sc_annotate")
sc_pseudobulk_dir = os.path.join(sc_dir, "sc_pseudobulk")
sc_dge_dir        = os.path.join(sc_dir, "sc_dge")
sc_pathway_dir    = os.path.join(sc_dir, "sc_pathway")

sc_markers = config.get("sc_markers", "resources/markers_mouse_brain.csv")
if not os.path.isabs(sc_markers):
    sc_markers = os.path.join(config.get("pipeline_dir"), sc_markers)
sc_threads        = int(config.get("sc_threads", 16))
sc_mem            = int(config.get("sc_mem", 200))          # GB, sc_cluster / sc_annotate
sc_plot_cols      = _csv("sc_plot_cols", "Condition 1,Sample sex,Batch")
sc_group_cols     = _csv("sc_group_cols", "Condition 1,Sample sex")
sc_celltypes      = _csv("sc_dge_celltypes")
sc_comparisons    = dict(config.get("sc_dge_comparisons", {}) or {})
for _c, _ref in sc_comparisons.items():
    if not _c.endswith("_vs_" + str(_ref)):
        raise ValueError(f"sc_dge_comparisons: '{_c}' must be NUMERATOR_vs_{_ref}")

wildcard_constraints:
    sample="|".join(re.escape(s) for s in samples),
    unit="|".join(re.escape(u) for u in units),
    celltype="[A-Za-z0-9_-]+",
    comparison="[A-Za-z0-9_.-]+"


# -----------------------------
# Rule all
# -----------------------------
rule all:
    input:
        # FastQC raw (one per FASTQ pair)
        expand(os.path.join(fastqc_dir, "{unit}.r_1_fastqc.html"), unit=units.keys()),
        expand(os.path.join(fastqc_dir, "{unit}.r_1_fastqc.zip"), unit=units.keys()),
        expand(os.path.join(fastqc_dir, "{unit}.r_2_fastqc.html"), unit=units.keys()),
        expand(os.path.join(fastqc_dir, "{unit}.r_2_fastqc.zip"), unit=units.keys()),

        # MultiQC raw
        os.path.join(multiqc_dir, "multiqc_report.html"),

        # QC check
        os.path.join(qc_check_dir, "qc_metrics.csv"),
        os.path.join(qc_check_dir, "qc_metrics_per_lane.csv"),
        os.path.join(qc_check_dir, "barplot_qc_run_lane.png"),
        [os.path.join(qc_check_dir, "stripplot_qc_groups.png")] if qc_check_cat else [],

        # Cell Ranger
        expand(os.path.join(cellranger_dir, "{sample}", "outs", "web_summary.html"), sample=samples.keys()),
        expand(os.path.join(cellranger_dir, "{sample}", "outs", "filtered_feature_bc_matrix.h5"), sample=samples.keys()),

        # Cell Ranger check
        os.path.join(cellranger_check_dir, "cellranger_metrics.csv"),
        os.path.join(cellranger_check_dir, "barplot_cellranger_metrics.png"),
        os.path.join(cellranger_check_dir, "barplot_cellranger_mapping.png"),

        # Pseudobulk counts
        os.path.join(merged_counts_dir, "merged_counts.csv"),

        # Sample check
        os.path.join(sample_check_dir, "heatmap_sex.png"),

        # Downstream single-cell analysis
        [os.path.join(sc_qc_check_dir, "qc_summary.csv"),
         os.path.join(sc_cluster_dir, "adata_clustered.h5ad"),
         os.path.join(sc_annotate_dir, "cells_annotated.csv.gz"),
         os.path.join(sc_dge_dir, "dge_summary.csv")] if sc_downstream else [],
        expand(os.path.join(sc_pathway_dir, "{celltype}", "{comparison}", "GO_GSEA.csv"),
               celltype=sc_celltypes, comparison=sc_comparisons) if sc_downstream else [],


# -----------------------------
# FastQC raw
# -----------------------------
rule fastqc:
    input:
        r1=lambda wc: units[wc.unit]["r1"],
        r2=lambda wc: units[wc.unit]["r2"]
    output:
        r1_html=os.path.join(fastqc_dir, "{unit}.r_1_fastqc.html"),
        r1_zip=os.path.join(fastqc_dir, "{unit}.r_1_fastqc.zip"),
        r2_html=os.path.join(fastqc_dir, "{unit}.r_2_fastqc.html"),
        r2_zip=os.path.join(fastqc_dir, "{unit}.r_2_fastqc.zip")
    threads:
        fastqc_threads,
    conda:
        "pipeline-bulkrnaseq-sm_env"
    shell:
        """
        mkdir -p {fastqc_dir}
        bash {SCRIPTS_DIR}/run_fastqc.sh {input.r1} {input.r2} {fastqc_dir} {threads}

        # For R1
        r1base=$(basename "{input.r1}" | sed -E 's/\\.(fq|fastq)(\\.gz)?$//')
        mv "{fastqc_dir}/${{r1base}}_fastqc.html" "{output.r1_html}"
        mv "{fastqc_dir}/${{r1base}}_fastqc.zip"  "{output.r1_zip}"

        # For R2
        r2base=$(basename "{input.r2}" | sed -E 's/\\.(fq|fastq)(\\.gz)?$//')
        mv "{fastqc_dir}/${{r2base}}_fastqc.html" "{output.r2_html}"
        mv "{fastqc_dir}/${{r2base}}_fastqc.zip"  "{output.r2_zip}"
        """

# -----------------------------
# MultiQC raw
# -----------------------------
rule multiqc:
    input:
        expand(os.path.join(fastqc_dir, "{unit}.r_1_fastqc.zip"), unit=units.keys()) +
        expand(os.path.join(fastqc_dir, "{unit}.r_2_fastqc.zip"), unit=units.keys())
    output:
        os.path.join(multiqc_dir, "multiqc_report.html")
    conda:
        "pipeline-bulkrnaseq-sm_env"
    shell:
        """
        mkdir -p {multiqc_dir}
        bash {SCRIPTS_DIR}/run_multiqc.sh {fastqc_dir} {multiqc_dir} {multiqc_mode}
        """

# -----------------------------
# QC check
# -----------------------------
# FastQC metrics per sample (lanes combined): a per-sample figure by sequencing run + lane and
# a comparison of all metrics (R1 and R2 side by side) across the qc_check_cat groups.
rule qc_check:
    input:
        multiqc_report=os.path.join(multiqc_dir, "multiqc_report.html"),
        filepaths=filepaths,
        metadata=metadata
    output:
        table=os.path.join(qc_check_dir, "qc_metrics.csv"),
        lane_table=os.path.join(qc_check_dir, "qc_metrics_per_lane.csv"),
        run_lane_plot=os.path.join(qc_check_dir, "barplot_qc_run_lane.png"),
        groups_plot=[os.path.join(qc_check_dir, "stripplot_qc_groups.png")] if qc_check_cat else []
    params:
        qc_check_cat=",".join(qc_check_cat)
    conda:
        "dge_env"
    shell:
        """
        mkdir -p {qc_check_dir}

        Rscript {SCRIPTS_DIR}/plot_barplot_qc.R \
            --multiqc_dir "{multiqc_dir}/multiqc_report_data" \
            --filepaths "{input.filepaths}" \
            --metadata "{input.metadata}" \
            --group_cols "{params.qc_check_cat}" \
            --outdir "{qc_check_dir}"
        """

# -----------------------------
# Cell Ranger count
# -----------------------------
# Waits on the raw MultiQC report (QC gate, as trim_reads in pipeline-bulkrnaseq-sm).
# Cell Ranger is not on conda: cellranger_bin is the full path to the 10x install.
# All FASTQ pairs of a sample are passed together (as lanes L001, L002, ...).
rule cellranger_count:
    input:
        r1=lambda wc: [p["r1"] for p in samples[wc.sample]],
        r2=lambda wc: [p["r2"] for p in samples[wc.sample]],
        multiqc_report=os.path.join(multiqc_dir, "multiqc_report.html")
    output:
        web=os.path.join(cellranger_dir, "{sample}", "outs", "web_summary.html"),
        metrics=os.path.join(cellranger_dir, "{sample}", "outs", "metrics_summary.csv"),
        filtered=os.path.join(cellranger_dir, "{sample}", "outs", "filtered_feature_bc_matrix.h5"),
        raw=os.path.join(cellranger_dir, "{sample}", "outs", "raw_feature_bc_matrix.h5")
    threads:
        cellranger_threads,
    resources:
        mem_mb=cellranger_mem * 1024
    params:
        r1=lambda wc: ",".join(p["r1"] for p in samples[wc.sample]),
        r2=lambda wc: ",".join(p["r2"] for p in samples[wc.sample])
    shell:
        """
        mkdir -p {cellranger_dir}
        bash {SCRIPTS_DIR}/run_cellranger.sh {cellranger_bin} {wildcards.sample} "{params.r1}" "{params.r2}" {cellranger_ref} {cellranger_dir} {threads} {cellranger_mem} {cellranger_chemistry} {cellranger_create_bam} {cellranger_introns} "{cellranger_expect_cells}"
        """

# -----------------------------
# Cell Ranger check
# -----------------------------
rule cellranger_check:
    input:
        metrics=expand(os.path.join(cellranger_dir, "{sample}", "outs", "metrics_summary.csv"), sample=samples.keys()),
        metadata=metadata
    output:
        table=os.path.join(cellranger_check_dir, "cellranger_metrics.csv"),
        metrics_plot=os.path.join(cellranger_check_dir, "barplot_cellranger_metrics.png"),
        mapping_plot=os.path.join(cellranger_check_dir, "barplot_cellranger_mapping.png")
    params:
        cellranger_check_cat=cellranger_check_cat,
        expected_cells=expected_cells
    conda:
        "pipeline-bulkrnaseq-sm_env"
    shell:
        """
        mkdir -p {cellranger_check_dir}

        python3 {SCRIPTS_DIR}/plot_cellranger_metrics.py \
            --metrics {input.metrics} \
            --metadata "{input.metadata}" \
            --group_col "{params.cellranger_check_cat}" \
            --expected_cells "{params.expected_cells}" \
            --outdir "{cellranger_check_dir}"
        """


# -----------------------------
# Pseudobulk counts
# -----------------------------
# Sum of UMIs over all cells per sample, same layout as merged_counts.csv in pipeline-bulkrnaseq-sm.
rule merge_counts:
    input:
        h5=expand(os.path.join(cellranger_dir, "{sample}", "outs", "filtered_feature_bc_matrix.h5"), sample=samples.keys())
    output:
        merged=os.path.join(merged_counts_dir, "merged_counts.csv"),
        ncells=os.path.join(merged_counts_dir, "n_cells.csv")
    conda:
        "pipeline-bulkrnaseq-sm_env"
    shell:
        """
        mkdir -p {merged_counts_dir}

        python3 {SCRIPTS_DIR}/make_pseudobulk.py \
            --h5 {input.h5} \
            --out "{output.merged}"
        """

# -----------------------------
# Sample QC Check (sex)
# -----------------------------
rule sex_check:
    input:
        counts=os.path.join(merged_counts_dir, "merged_counts.csv"),
        metadata=metadata,
        annotation=annotation
    output:
        sex_plot=os.path.join(sample_check_dir, "heatmap_sex.png")
    conda:
        "dge_env"
    shell:
        """
        mkdir -p {sample_check_dir}

        Rscript {SCRIPTS_DIR}/plot_sex.R \
            --counts "{input.counts}" \
            --metadata "{input.metadata}" \
            --annotation "{input.annotation}" \
            --out "{output.sex_plot}"
        """


# =============================================================================
# Downstream single-cell analysis
# =============================================================================
# sc_qc -> sc_qc_check (gate) -> sc_cluster -> sc_annotate -> sc_dge -> sc_pathway, sc_dge_summary

# -----------------------------
# Cell QC per sample (QC metrics + Scrublet doublets)
# -----------------------------
rule sc_qc:
    input:
        h5=os.path.join(cellranger_dir, "{sample}", "outs", "filtered_feature_bc_matrix.h5")
    output:
        cells=os.path.join(sc_qc_dir, "{sample}_cells.csv.gz")
    params:
        min_genes=config.get("sc_min_genes", 300),
        max_genes=config.get("sc_max_genes", 10000),
        min_counts=config.get("sc_min_counts", 500),
        max_pct_mito=config.get("sc_max_pct_mito", 5),
        doublets=str(config.get("sc_doublets", True)).lower(),
        doublet_rate=config.get("sc_doublet_rate_per_1k", 0.004)
    resources:
        mem_mb=16000
    conda:
        "pipeline-scrnaseq-sm_env"
    shell:
        """
        python3 {SCRIPTS_DIR}/sc_qc.py \
            --h5 "{input.h5}" --sample {wildcards.sample} --out "{output.cells}" \
            --min_genes {params.min_genes} --max_genes {params.max_genes} \
            --min_counts {params.min_counts} --max_pct_mito {params.max_pct_mito} \
            --doublets {params.doublets} --doublet_rate_per_1k {params.doublet_rate}
        """

# -----------------------------
# Cell QC across samples (gate: clustering waits for these plots)
# -----------------------------
rule sc_qc_check:
    input:
        cells=expand(os.path.join(sc_qc_dir, "{sample}_cells.csv.gz"), sample=samples.keys()),
        metadata=metadata
    output:
        summary=os.path.join(sc_qc_check_dir, "qc_summary.csv"),
        violin=os.path.join(sc_qc_check_dir, "violin_qc_per_sample.png"),
        barplot=os.path.join(sc_qc_check_dir, "barplot_cell_status.png")
    params:
        group_cols=",".join(sc_group_cols),
        min_genes=config.get("sc_min_genes", 300),
        max_genes=config.get("sc_max_genes", 10000),
        min_counts=config.get("sc_min_counts", 500),
        max_pct_mito=config.get("sc_max_pct_mito", 5)
    conda:
        "pipeline-scrnaseq-sm_env"
    shell:
        """
        python3 {SCRIPTS_DIR}/sc_qc_check.py \
            --qc {input.cells} --metadata "{input.metadata}" --group_cols "{params.group_cols}" \
            --min_genes {params.min_genes} --max_genes {params.max_genes} \
            --min_counts {params.min_counts} --max_pct_mito {params.max_pct_mito} \
            --outdir "{sc_qc_check_dir}"
        """

# -----------------------------
# Normalise, integrate (Harmony), UMAP, Leiden, cluster markers
# -----------------------------
rule sc_cluster:
    input:
        h5=expand(os.path.join(cellranger_dir, "{sample}", "outs", "filtered_feature_bc_matrix.h5"), sample=samples.keys()),
        cells=expand(os.path.join(sc_qc_dir, "{sample}_cells.csv.gz"), sample=samples.keys()),
        qc_gate=os.path.join(sc_qc_check_dir, "qc_summary.csv"),
        metadata=metadata
    output:
        h5ad=os.path.join(sc_cluster_dir, "adata_clustered.h5ad"),
        markers=os.path.join(sc_cluster_dir, "cluster_markers.csv"),
        umap=os.path.join(sc_cluster_dir, "umap_clusters.png"),
        umap_meta=os.path.join(sc_cluster_dir, "umap_metadata.png"),
        pca=os.path.join(sc_cluster_dir, "pca_pseudobulk_samples.png")
    params:
        n_hvg=config.get("sc_n_hvg", 3000),
        n_pcs=config.get("sc_n_pcs", 30),
        batch_key=config.get("sc_batch_key", "sample"),
        n_neighbors=config.get("sc_n_neighbors", 15),
        resolutions=",".join(_csv("sc_resolutions", "0.2,0.5,1.0")),
        cluster_resolution=config.get("sc_cluster_resolution", 0.5),
        plot_cols=",".join(sc_plot_cols),
        split_cols=",".join(sc_group_cols),
        subsample=config.get("sc_subsample_cells", "") or ""
    threads:
        sc_threads
    resources:
        mem_mb=sc_mem * 1024
    conda:
        "pipeline-scrnaseq-sm_env"
    shell:
        """
        python3 {SCRIPTS_DIR}/sc_cluster.py \
            --h5 {input.h5} --qc {input.cells} --metadata "{input.metadata}" --outdir "{sc_cluster_dir}" \
            --n_hvg {params.n_hvg} --n_pcs {params.n_pcs} --batch_key "{params.batch_key}" \
            --n_neighbors {params.n_neighbors} --resolutions "{params.resolutions}" \
            --cluster_resolution {params.cluster_resolution} --plot_cols "{params.plot_cols}" \
            --split_cols "{params.split_cols}" --subsample "{params.subsample}" --threads {threads}
        """

# -----------------------------
# Cell types, subclustering, composition, pseudobulk per cell type
# -----------------------------
rule sc_annotate:
    input:
        h5ad=os.path.join(sc_cluster_dir, "adata_clustered.h5ad"),
        markers=sc_markers,
        metadata=metadata
    output:
        cells=os.path.join(sc_annotate_dir, "cells_annotated.csv.gz"),
        umap=os.path.join(sc_annotate_dir, "umap_celltype.png"),
        composition=os.path.join(sc_annotate_dir, "composition.csv"),
        ncells=os.path.join(sc_pseudobulk_dir, "n_cells_per_sample.csv"),
        counts=expand(os.path.join(sc_pseudobulk_dir, "{celltype}", "counts.csv"), celltype=sc_celltypes),
        pb_meta=expand(os.path.join(sc_pseudobulk_dir, "{celltype}", "metadata.csv"), celltype=sc_celltypes),
        status=expand(os.path.join(sc_pseudobulk_dir, "{celltype}", "status_{comparison}.txt"),
                      celltype=sc_celltypes, comparison=sc_comparisons)
    params:
        plot_cols=",".join(sc_plot_cols),
        group_cols=",".join(sc_group_cols),
        subcluster_types=",".join(_csv("sc_subcluster_types")),
        subcluster_resolution=config.get("sc_subcluster_resolution", 0.5),
        celltypes=",".join(sc_celltypes),
        comparisons=",".join(sc_comparisons),
        min_cells=config.get("sc_min_cells_pseudobulk", 10)
    threads:
        sc_threads
    resources:
        mem_mb=sc_mem * 1024
    conda:
        "pipeline-scrnaseq-sm_env"
    shell:
        """
        python3 {SCRIPTS_DIR}/sc_annotate.py \
            --h5ad "{input.h5ad}" --markers "{input.markers}" --metadata "{input.metadata}" \
            --outdir "{sc_annotate_dir}" --pseudobulk_dir "{sc_pseudobulk_dir}" \
            --plot_cols "{params.plot_cols}" --group_cols "{params.group_cols}" \
            --subcluster_types "{params.subcluster_types}" --subcluster_resolution {params.subcluster_resolution} \
            --celltypes "{params.celltypes}" --comparisons "{params.comparisons}" \
            --min_cells {params.min_cells} --threads {threads}
        """

# -----------------------------
# Pseudobulk DESeq2 per cell type and comparison (bulk run_deseq_dge.R, design ~ Group)
# -----------------------------
# Skipped (empty dge_results.csv) when a group has < 2 samples with enough cells of the type.
rule sc_dge:
    input:
        counts=os.path.join(sc_pseudobulk_dir, "{celltype}", "counts.csv"),
        metadata=os.path.join(sc_pseudobulk_dir, "{celltype}", "metadata.csv"),
        status=os.path.join(sc_pseudobulk_dir, "{celltype}", "status_{comparison}.txt"),
        annotation=annotation
    output:
        results=os.path.join(sc_dge_dir, "{celltype}", "{comparison}", "dge_results.csv")
    params:
        reference=lambda wc: sc_comparisons[wc.comparison],
        lfc_threshold=config.get("lfc_threshold", 1.0),
        padj_threshold=config.get("padj_threshold", 0.05)
    conda:
        "dge_env"
    shell:
        """
        mkdir -p "$(dirname {output.results})"
        if [ "$(cat {input.status})" = "ok" ]; then
            Rscript {SCRIPTS_DIR}/run_deseq_dge.R \
                --counts "{input.counts}" --metadata "{input.metadata}" --annotation "{input.annotation}" \
                --out "{output.results}" --comparison "{wildcards.comparison}" --reference "{params.reference}" \
                --dge_cat "Group" --filt TRUE \
                --lfc_threshold {params.lfc_threshold} --padj_threshold {params.padj_threshold}
        else
            echo "[INFO] {wildcards.celltype} {wildcards.comparison}: $(cat {input.status})"
            : > "{output.results}"
        fi
        """

# -----------------------------
# GSEA per cell type and comparison (bulk run_gsea_pathway.R)
# -----------------------------
rule sc_pathway:
    input:
        dge=os.path.join(sc_dge_dir, "{celltype}", "{comparison}", "dge_results.csv"),
        annotation=annotation
    output:
        go=os.path.join(sc_pathway_dir, "{celltype}", "{comparison}", "GO_GSEA.csv"),
        kegg=os.path.join(sc_pathway_dir, "{celltype}", "{comparison}", "KEGG_GSEA.csv"),
        reactome=os.path.join(sc_pathway_dir, "{celltype}", "{comparison}", "REACTOME_GSEA.csv")
    params:
        outdir=lambda wc: os.path.join(sc_pathway_dir, wc.celltype, wc.comparison),
        organism=config.get("pathway_organism", "mouse"),
        databases=",".join(_csv("pathway_databases", "GO,KEGG,REACTOME")),
        min_genes=config.get("pathway_min_genes", 10),
        max_genes=config.get("pathway_max_genes", 500),
        gsea_fdr=config.get("pathway_gsea_fdr", 0.05)
    conda:
        "dge_env"
    shell:
        """
        mkdir -p "{params.outdir}"
        if [ -s "{input.dge}" ]; then
            Rscript {SCRIPTS_DIR}/run_gsea_pathway.R \
                --dge "{input.dge}" --annotation "{input.annotation}" --outdir "{params.outdir}" \
                --comparison "{wildcards.comparison}" --organism "{params.organism}" \
                --databases "{params.databases}" --min_genes {params.min_genes} \
                --max_genes {params.max_genes} --gsea_fdr {params.gsea_fdr}
        fi
        for f in {output.go} {output.kegg} {output.reactome}; do [ -e "$f" ] || : > "$f"; done
        """

# -----------------------------
# DESeq2 summary across cell types
# -----------------------------
rule sc_dge_summary:
    input:
        dge=expand(os.path.join(sc_dge_dir, "{celltype}", "{comparison}", "dge_results.csv"),
                   celltype=sc_celltypes, comparison=sc_comparisons),
        ncells=os.path.join(sc_pseudobulk_dir, "n_cells_per_sample.csv")
    output:
        table=os.path.join(sc_dge_dir, "dge_summary.csv"),
        plot=os.path.join(sc_dge_dir, "barplot_dge_summary.png")
    params:
        celltypes=",".join(sc_celltypes),
        comparisons=",".join(sc_comparisons)
    conda:
        "pipeline-scrnaseq-sm_env"
    shell:
        """
        python3 {SCRIPTS_DIR}/sc_dge_summary.py \
            --dge_dir "{sc_dge_dir}" --pseudobulk_dir "{sc_pseudobulk_dir}" \
            --celltypes "{params.celltypes}" --comparisons "{params.comparisons}" --outdir "{sc_dge_dir}"
        """

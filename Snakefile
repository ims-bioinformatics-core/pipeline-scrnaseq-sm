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

        if sample in samples:
            raise ValueError(f"Sample '{sample}' appears more than once in {filepaths}. Use one row (one R1/R2 pair) per sample - merge lanes first.")

        fq1 = os.path.join(RUN_DIR, parts[1])
        fq2 = os.path.join(RUN_DIR, parts[2])
        samples[sample] = {"r1": fq1, "r2": fq2}

print("[INFO] Samples:", ", ".join(samples))

metadata = os.path.join(RUN_DIR, config["metadata"])
# -----------------------------
# Directories for outputs under RUN_DIR
# -----------------------------

DIR_FLAG = config["dir_flag"]

fastqc_dir        = os.path.join(RUN_DIR, DIR_FLAG, "fastq-raw")
multiqc_dir       = os.path.join(RUN_DIR, DIR_FLAG, "multiqc-raw")
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


# -----------------------------
# Rule all
# -----------------------------
rule all:
    input:
        # FastQC raw
        expand(os.path.join(fastqc_dir, "{sample}.r_1_fastqc.html"), sample=samples.keys()),
        expand(os.path.join(fastqc_dir, "{sample}.r_1_fastqc.zip"), sample=samples.keys()),
        expand(os.path.join(fastqc_dir, "{sample}.r_2_fastqc.html"), sample=samples.keys()),
        expand(os.path.join(fastqc_dir, "{sample}.r_2_fastqc.zip"), sample=samples.keys()),

        # MultiQC raw
        os.path.join(multiqc_dir, "multiqc_report.html"),

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


# -----------------------------
# FastQC raw
# -----------------------------
rule fastqc:
    input:
        r1=lambda wc: samples[wc.sample]["r1"],
        r2=lambda wc: samples[wc.sample]["r2"]
    output:
        r1_html=os.path.join(fastqc_dir, "{sample}.r_1_fastqc.html"),
        r1_zip=os.path.join(fastqc_dir, "{sample}.r_1_fastqc.zip"),
        r2_html=os.path.join(fastqc_dir, "{sample}.r_2_fastqc.html"),
        r2_zip=os.path.join(fastqc_dir, "{sample}.r_2_fastqc.zip")
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
        expand(os.path.join(fastqc_dir, "{sample}.r_1_fastqc.zip"), sample=samples.keys()) +
        expand(os.path.join(fastqc_dir, "{sample}.r_2_fastqc.zip"), sample=samples.keys())
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
# Cell Ranger count
# -----------------------------
# Waits on the raw MultiQC report (QC gate, as trim_reads in pipeline-bulkrnaseq-sm).
# Cell Ranger is not on conda: cellranger_bin is the full path to the 10x install.
rule cellranger_count:
    input:
        r1=lambda wc: samples[wc.sample]["r1"],
        r2=lambda wc: samples[wc.sample]["r2"],
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
    shell:
        """
        mkdir -p {cellranger_dir}
        bash {SCRIPTS_DIR}/run_cellranger.sh {cellranger_bin} {wildcards.sample} {input.r1} {input.r2} {cellranger_ref} {cellranger_dir} {threads} {cellranger_mem} {cellranger_chemistry} {cellranger_create_bam} {cellranger_introns} "{cellranger_expect_cells}"
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

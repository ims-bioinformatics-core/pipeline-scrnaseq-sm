# pipeline-scrnaseq-sm

Single-cell RNA-seq (10x Genomics) pipeline in Snakemake, built the same way as `pipeline-bulkrnaseq-sm`
(same `run.sh` / `config.yml` / `inputs/` / `scripts/` layout, same `filepaths.csv` and `metadata.csv`).

## Workflow

| # | Rule | Tool | Env | Produces |
|---|------|------|-----|----------|
| 1 | `fastqc` | FastQC | `pipeline-bulkrnaseq-sm_env` | `fastq-raw/{sample}.r_{1,2}_fastqc.{html,zip}` |
| 2 | `multiqc` | MultiQC | `pipeline-bulkrnaseq-sm_env` | `multiqc-raw/multiqc_report.html` |
| 3 | `cellranger_count` | Cell Ranger 10.1.0 | none (full path in `cellranger_bin`) | `cellranger/{sample}/outs/` |
| 4 | `cellranger_check` | Python (`plot_cellranger_metrics.py`) | `pipeline-bulkrnaseq-sm_env` | `cellranger_check/`: merged metrics table + per-sample QC barplots, grouped by `cellranger_check_cat` |
| 5 | `merge_counts` | Python (`make_pseudobulk.py`) | `pipeline-bulkrnaseq-sm_env` | `merged-counts/merged_counts.csv`: pseudobulk (UMIs summed over cells) gene x sample, same layout as bulk; `n_cells.csv` |
| 6 | `sex_check` | R (`plot_sex.R`, from bulk) | `dge_env` | `sample_check/heatmap_sex.png`: Xist vs Y-gene expression per sample, annotated with metadata sex |

As in the bulk pipeline, the raw MultiQC report is a gate: Cell Ranger starts only after it exists.

## Requirements

```bash
# FastQC + MultiQC run in the existing pipeline-bulkrnaseq-sm_env (referenced by name in the Snakefile)
```

Cell Ranger is not on conda: `config.yml` points to `/rds/project/rds-O11U8YqSuCk/bioinformatics/tools/cellranger-10.1.0/cellranger`.
Reference: Ensembl GRCm39 release 113 (`genome/index/cellranger/ensembl/mus_musculus/release-113/mus_musculus_GRCm39_113`, built with cellranger mkref 10.1.0).

## Inputs

- `inputs/filepaths.csv` — `sample,fastq_1,fastq_2`, **one row per sample** (paths relative to `run_dir`).
  R1 = barcode + UMI read, R2 = cDNA read. Sample names: letters, digits, `-`, `_` only.
  Generate with `bash prerun_filepaths.sh` (`./raw_data` → `./inputs/filepaths.csv`).
- `inputs/metadata.csv` — same format as bulk (`Barcode` column = sample names). Not used by these three steps; kept for downstream analysis.

`raw_data` is a symlink to `../pipeline-scrnaseq-sm_old/raw_data`.

## Running

```bash
sbatch submit_slurm.sh      # 24 cores / 78G / 12h on icelake -> FastQC 12 samples x 2 threads, then Cell Ranger 1 sample at a time (24 cores, 64 GB)
# or interactively
bash run.sh
```

Thread / memory settings are at the top of `run.sh` (`CELLRANGER_THREADS`, `CELLRANGER_MEM`, `FASTQC_THREADS`);
cores come from `SLURM_CPUS_PER_TASK`. If the job hits the time limit, resubmit: finished samples are kept and
unfinished Cell Ranger runs restart from scratch (`--rerun-incomplete`). If Snakemake reports a lock after a killed job,
run once with `--unlock`.

```bash
# Dry run
snakemake -s Snakefile -n --cores 96 --configfile config.yml \
  --config cellranger_threads=16 cellranger_mem=64 fastqc_threads=4 pipeline_dir=$PWD
# QC only
... --until multiqc
# Skip samples
... --config exclude="SITTA5,SITTB5"
```

## Config (`config.yml`)

| Key | Description |
|-----|-------------|
| `dir_flag` | Output folder under `run_dir` (`results_raw`) |
| `run_dir`, `filepaths`, `metadata` | As in the bulk pipeline |
| `multiqc_mode` | `interactive` or `flat` |
| `exclude` | Optional comma-separated samples to skip |
| `cellranger_bin` | Full path to `cellranger` |
| `cellranger_ref` | Cell Ranger reference directory |
| `cellranger_chemistry` | `auto` or e.g. `SC3Pv3`, `SC3Pv4`, `SC5P-R2` |
| `cellranger_create_bam` | `true` / `false` |
| `cellranger_include_introns` | `true` (default) / `false` |
| `cellranger_expect_cells` | Optional integer, empty = auto |
| `cellranger_check_cat` | Metadata column to group the Cell Ranger QC plots by (e.g. `Sample sex`; empty = no grouping) |
| `annotation` | Gene annotation TSV (`ensembl_gene_id`, `gene_symbol`), same Ensembl release as `cellranger_ref`; used by `sex_check` |
| `expected_cells` | Plot only: dashed line at this value on "Estimated Number of Cells" (empty = no line). Does not change Cell Ranger (use `cellranger_expect_cells` for that) |

## Outputs

```
results_raw/
├── fastq-raw/                     # FastQC per sample
├── multiqc-raw/multiqc_report.html
├── cellranger/
│   ├── fastq/<sample>/            # symlinks named <sample>_S1_L001_R{1,2}_001.fastq.gz for Cell Ranger
│   └── <sample>/outs/             # web_summary.html, metrics_summary.csv, *_feature_bc_matrix(.h5), BAM, cloupe
├── merged-counts/
│   ├── merged_counts.csv          # pseudobulk counts: gene_id x sample (UMIs summed over all cells)
│   └── n_cells.csv                # cells per sample in the filtered matrix
├── sample_check/
│   └── heatmap_sex.png            # sex-marker heatmap (as in pipeline-bulkrnaseq-sm)
└── cellranger_check/
    ├── cellranger_metrics.csv           # one row per sample: metadata + all Cell Ranger metrics as numbers
    ├── barplot_cellranger_metrics.png   # 12 key metrics per sample, dashed lines = 10x guidelines
    └── barplot_cellranger_mapping.png   # exonic / intronic / intergenic / unmapped reads per sample
```

## Author

Maintained by mk2314, University of Cambridge

#!/bin/bash

# ==========================
# Pipeline directories
# ==========================
PIPELINE_DIR="/rds/project/rds-O11U8YqSuCk/bioinformatics/pipelines/pipeline-scrnaseq-sm"
CONFIG_FILE="$PIPELINE_DIR/config.yml"

# ==========================
# Resources
# ==========================
# Everything is sized from the actual allocation, so only submit_slurm.sh needs editing.
# nproc = CPUs really allocated (SLURM may add CPUs to cover --mem); memory from SLURM, in GB.
CORES=$(nproc)
if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
    TOTAL_MEM=$((SLURM_MEM_PER_NODE / 1024))
elif [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
    TOTAL_MEM=$((SLURM_MEM_PER_CPU * CORES / 1024))
else
    TOTAL_MEM=64                         # interactive default
fi

CELLRANGER_MAX_THREADS=24                # cores per Cell Ranger job (capped at CORES)
CELLRANGER_MAX_MEM=64                    # GB per Cell Ranger job (capped at TOTAL_MEM - 8)
FASTQC_THREADS=$(( CORES < 2 ? CORES : 2 ))   # FastQC uses 1 thread per file -> 2 per sample (R1 + R2) is the max useful

CELLRANGER_THREADS=$(( CORES < CELLRANGER_MAX_THREADS ? CORES : CELLRANGER_MAX_THREADS ))
CELLRANGER_MEM=$(( TOTAL_MEM - 8 < CELLRANGER_MAX_MEM ? TOTAL_MEM - 8 : CELLRANGER_MAX_MEM ))
if [ "$CELLRANGER_MEM" -lt 32 ]; then
    echo "❌ Only ${TOTAL_MEM} GB allocated; Cell Ranger needs >= 40 GB (mouse/human). Increase --mem."
    exit 1
fi
PARALLEL_CR=$(( CORES / CELLRANGER_THREADS < TOTAL_MEM / CELLRANGER_MEM ? CORES / CELLRANGER_THREADS : TOTAL_MEM / CELLRANGER_MEM ))

echo "🖥️  Allocation: ${CORES} cores, ${TOTAL_MEM} GB"
echo "    Cell Ranger: ${CELLRANGER_THREADS} cores / ${CELLRANGER_MEM} GB per sample, ${PARALLEL_CR} sample(s) at a time"
echo "    FastQC:      ${FASTQC_THREADS} threads"

# ==========================
# Run Snakemake
# ==========================
# Clear a stale lock left by a job killed at the time limit (safe: one pipeline job at a time).
CONFIG_ARGS=(cellranger_threads=$CELLRANGER_THREADS cellranger_mem=$CELLRANGER_MEM fastqc_threads=$FASTQC_THREADS pipeline_dir="$PIPELINE_DIR")
snakemake -s "$PIPELINE_DIR/Snakefile" --unlock --configfile "$CONFIG_FILE" --config "${CONFIG_ARGS[@]}" > /dev/null 2>&1

start=$(date +%s)
snakemake -s "$PIPELINE_DIR/Snakefile" \
          --cores "$CORES" \
          --resources mem_mb=$((TOTAL_MEM * 1024)) \
          --use-conda \
          --rerun-incomplete \
          --configfile "$CONFIG_FILE" \
          --config "${CONFIG_ARGS[@]}"
status=$?

end=$(date +%s)
runtime=$((end - start))
runtime_minutes=$(echo "scale=2; $runtime/60" | bc)
if [ $status -eq 0 ]; then
    echo "Pipeline finished in $runtime seconds ($runtime_minutes minutes)"
else
    echo "❌ Pipeline FAILED (exit $status) after $runtime seconds ($runtime_minutes minutes)"
fi
exit $status

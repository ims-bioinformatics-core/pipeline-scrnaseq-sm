#!/bin/bash
# Usage: run_cellranger.sh <cellranger_bin> <sample> <read1> <read2> <transcriptome> <out_dir>
#                          <threads> <mem_gb> <chemistry> <create_bam> <include_introns> [expect_cells]
set -euo pipefail


# --- Inputs ---
CELLRANGER="$1"
SAMPLE="$2"
R1=$(realpath "$3")
R2=$(realpath "$4")
TRANSCRIPTOME=$(realpath "$5")
OUT_DIR=$(realpath "$6")
THREADS="$7"
MEM_GB="$8"
CHEMISTRY="$9"
CREATE_BAM="${10}"
INCLUDE_INTRONS="${11}"
EXPECT_CELLS="${12:-}"

# --- Prepare FASTQs ---
# Cell Ranger needs bcl2fastq-style names: {sample}_S1_L001_R{1,2}_001.fastq.gz
FASTQ_DIR="$OUT_DIR/fastq/$SAMPLE"
mkdir -p "$FASTQ_DIR"
ln -sf "$R1" "$FASTQ_DIR/${SAMPLE}_S1_L001_R1_001.fastq.gz"
ln -sf "$R2" "$FASTQ_DIR/${SAMPLE}_S1_L001_R2_001.fastq.gz"

# --- Prepare output directory ---
# Snakemake pre-creates <sample>/outs/ for the declared outputs, and Cell Ranger refuses
# to write into a directory that is not one of its own pipestances -> start clean.
cd "$OUT_DIR"
rm -rf "$OUT_DIR/$SAMPLE"

OPTIONAL_ARGS=()
[[ -n "$EXPECT_CELLS" ]] && OPTIONAL_ARGS+=(--expect-cells="$EXPECT_CELLS")

echo "🚀 Cell Ranger count:"
echo "    Cell Ranger: $("$CELLRANGER" --version 2>&1 | head -1) ($CELLRANGER)"
echo "    Sample: $SAMPLE"
echo "    R1: $R1"
echo "    R2: $R2"
echo "    Transcriptome: $TRANSCRIPTOME"
echo "    Output dir: $OUT_DIR/$SAMPLE"
echo "    Cores / mem: $THREADS / ${MEM_GB} GB"
echo "    Chemistry: $CHEMISTRY"
echo "    Create BAM: $CREATE_BAM"
echo "    Include introns: $INCLUDE_INTRONS"
echo "    Expect cells: ${EXPECT_CELLS:-auto}"

# --- Run Cell Ranger ---
"$CELLRANGER" count \
    --id="$SAMPLE" \
    --sample="$SAMPLE" \
    --fastqs="$FASTQ_DIR" \
    --transcriptome="$TRANSCRIPTOME" \
    --chemistry="$CHEMISTRY" \
    --create-bam="$CREATE_BAM" \
    --include-introns="$INCLUDE_INTRONS" \
    --localcores="$THREADS" \
    --localmem="$MEM_GB" \
    --disable-ui \
    ${OPTIONAL_ARGS[@]+"${OPTIONAL_ARGS[@]}"}

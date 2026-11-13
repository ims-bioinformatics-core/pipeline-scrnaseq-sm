#!/bin/bash
# Usage: run_cellranger.sh <cellranger_bin> <sample> <read1[,read1,...]> <read2[,read2,...]> <transcriptome> <out_dir>
#                          <threads> <mem_gb> <chemistry> <create_bam> <include_introns> [expect_cells]
# Several lanes / runs of one sample: comma-separated R1 and R2 lists in the same order.
set -euo pipefail


# --- Inputs ---
CELLRANGER="$1"
SAMPLE="$2"
IFS=',' read -r -a R1_LIST <<< "$3"
IFS=',' read -r -a R2_LIST <<< "$4"
TRANSCRIPTOME=$(realpath "$5")
OUT_DIR=$(realpath "$6")
THREADS="$7"
MEM_GB="$8"
CHEMISTRY="$9"
CREATE_BAM="${10}"
INCLUDE_INTRONS="${11}"
EXPECT_CELLS="${12:-}"

if [ "${#R1_LIST[@]}" -ne "${#R2_LIST[@]}" ]; then
    echo "ERROR: ${#R1_LIST[@]} R1 files but ${#R2_LIST[@]} R2 files for $SAMPLE"
    exit 1
fi

# --- Prepare FASTQs ---
# Cell Ranger needs bcl2fastq-style names: {sample}_S1_L00N_R{1,2}_001.fastq.gz (one "lane" per FASTQ pair)
FASTQ_DIR="$OUT_DIR/fastq/$SAMPLE"
rm -rf "$FASTQ_DIR"
mkdir -p "$FASTQ_DIR"
for i in "${!R1_LIST[@]}"; do
    LANE=$(printf "L%03d" $((i + 1)))
    ln -s "$(realpath "${R1_LIST[$i]}")" "$FASTQ_DIR/${SAMPLE}_S1_${LANE}_R1_001.fastq.gz"
    ln -s "$(realpath "${R2_LIST[$i]}")" "$FASTQ_DIR/${SAMPLE}_S1_${LANE}_R2_001.fastq.gz"
done

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
echo "    FASTQ pairs: ${#R1_LIST[@]}"
for i in "${!R1_LIST[@]}"; do
    echo "      L$(printf "%03d" $((i + 1))): ${R1_LIST[$i]}"
    echo "            ${R2_LIST[$i]}"
done
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

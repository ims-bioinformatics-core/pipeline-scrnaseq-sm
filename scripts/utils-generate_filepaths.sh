#!/bin/bash
# Build inputs/filepaths.csv (sample,fastq_1,fastq_2) from a directory of 10x FASTQs.
# One row per R1/R2 pair; a sample sequenced on several lanes / flowcells gets several rows,
# which the pipeline passes together to a single `cellranger count` run.
#
# Recognised naming:
#   CRUK Genomics :  SLX-12345.SIGAA1.HXXXXXXXX.s_1.r_1.fq.gz   -> sample = SIGAA1 (2nd dot field)
#   bcl2fastq     :  MySample_S1_L001_R1_001.fastq.gz           -> sample = MySample
#   generic       :  barcode_sample_1.fq.gz                     -> sample = barcode_sample
# Index reads (I1/I2, r_3/r_4) are ignored - Cell Ranger count does not need them.

if [ $# -ne 3 ]; then
    echo "Usage: $0 <input_directory> <output_directory> <output_csv_filename>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
CSV_FILENAME="$3"

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

RUN_DIR=$(pwd)
mkdir -p "$OUTPUT_DIR"
OUTPUT_FILE="$OUTPUT_DIR/$CSV_FILENAME"
TMP_FILE="${OUTPUT_FILE}.tmp"
> "$TMP_FILE"

echo "Searching FASTQ files in: $INPUT_DIR"
echo "Saving relative paths from: $RUN_DIR"

find "$INPUT_DIR" -type f \( -name "*.r_1.fq.gz" -o -name "*_R1_001.fastq.gz" -o -name "*_1.fq.gz" \) | sort | while read -r R1; do

    BASENAME=$(basename "$R1")

    if [[ "$BASENAME" == *.r_1.fq.gz ]]; then
        R2="${R1%.r_1.fq.gz}.r_2.fq.gz"
        SAMPLE=$(echo "$BASENAME" | cut -d'.' -f2)

    elif [[ "$BASENAME" =~ ^(.+)_S[0-9]+_L[0-9]{3}_R1_001\.fastq\.gz$ ]]; then
        R2="${R1%_R1_001.fastq.gz}_R2_001.fastq.gz"
        SAMPLE="${BASH_REMATCH[1]}"

    elif [[ "$BASENAME" =~ _1\.fq\.gz$ ]]; then
        R2="${R1%_1.fq.gz}_2.fq.gz"
        SAMPLE=$(echo "$BASENAME" | awk -F'_' '{print $1"_"$2}')

    else
        continue
    fi

    if [ -f "$R2" ]; then
        R1_PATH=$(realpath --relative-to="$RUN_DIR" "$R1")
        R2_PATH=$(realpath --relative-to="$RUN_DIR" "$R2")
        echo "$SAMPLE,$R1_PATH,$R2_PATH" >> "$TMP_FILE"
    else
        echo "Warning: Missing R2 file for $R1" >&2
    fi
done

{
    echo "Sample,fastq_1,fastq_2"
    sort -t',' -k1,1 -k2,2 "$TMP_FILE"
} > "$OUTPUT_FILE"

rm -f "$TMP_FILE"

echo "CSV file generated: $OUTPUT_FILE"
echo "Samples and number of FASTQ pairs:"
tail -n +2 "$OUTPUT_FILE" | cut -d',' -f1 | sort | uniq -c

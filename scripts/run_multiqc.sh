#!/bin/bash
# Usage: ./run_multiqc.sh <fastqc_dir> <multiqc_dir> <interactive|flat>
set -euo pipefail

FASTQC_DIR=$1
MULTIQC_DIR=$2
MULTIQC_MODE="$3"

# Choose mode
if [[ "$MULTIQC_MODE" == "interactive" ]]; then
    INTERACTIVE_FLAG="--interactive"
elif [[ "$MULTIQC_MODE" == "flat" ]]; then
    INTERACTIVE_FLAG="--flat"
else
    echo "ERROR: mode must be 'interactive' or 'flat'"
    exit 1
fi

mkdir -p "$MULTIQC_DIR"

multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR" -n "multiqc_report.html" -f $INTERACTIVE_FLAG

#!/bin/bash
# Usage: run_fastqc.sh read1 read2 outdir threads
set -euo pipefail

R1=$1
R2=$2
OUTDIR=$3
THREADS=$4

fastqc -o "$OUTDIR" -f fastq "$R1" "$R2" -t "$THREADS"

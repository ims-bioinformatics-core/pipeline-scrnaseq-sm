#!/bin/bash
# Split + downsample neuron_1k_v3 into 4 pseudo-samples (keep 20% of read pairs), both lanes in parallel
set -euo pipefail
T=/rds/project/rds-O11U8YqSuCk/bioinformatics/pipelines/pipeline-scrnaseq-sm/test_data
export PATH=/rds/project/rds-O11U8YqSuCk/bioinformatics/envs/miniconda3/envs/pipeline-bulkrnaseq-sm_env/bin:$PATH
IN=$T/download/neuron_1k_v3_fastqs
for l in 1 2; do
  python3 $T/split_by_barcode.py $IN/neuron_1k_v3_S1_L00${l}_R1_001.fastq.gz $IN/neuron_1k_v3_S1_L00${l}_R2_001.fastq.gz $T/fastq $l 4 0.2 &
done
wait
echo ALL_DONE

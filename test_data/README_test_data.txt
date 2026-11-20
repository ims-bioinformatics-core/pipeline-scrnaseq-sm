Test data for pipeline-scrnaseq-sm
==================================

Source (public, 10x Genomics):
  "1k Brain Cells from an E18 Mouse (v3 chemistry)", Cell Ranger 3.0.0 example dataset
  https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v3/neuron_1k_v3_fastqs.tar
  Mouse (matches the pipeline's GRCm39 / Ensembl 113 Cell Ranger reference), 3' v3, ~1,200 cells, 2 lanes.

Pseudo-samples:
  The library is split by cell barcode (first 16 bases of R1) into 4 pseudo-samples TEST1-TEST4
  with split_by_barcode.py, so each cell is in exactly one pseudo-sample. Each lane is kept as its own
  FASTQ pair (s_1, s_2), so every sample has 2 pairs (tests multi-lane handling).

  test_data/download/   original tar + extracted FASTQs (not in git)
  test_data/fastq/      SLX-00000.TEST<k>.NEURON1K.s_<lane>.r_{1,2}.fq.gz (not in git)

Metadata (inputs/metadata.csv) is FAKE: Condition 1 (treated/untreated), Condition 2, Batch,
Sample sex are made up so that every step runs (DGE needs 2 samples per group). The real
sample is a pool of E18 embryonic brains, so sex_check and DGE results have no biological meaning.

To rebuild from scratch:
  cd test_data/download && curl -LO <URL above> && tar -xf neuron_1k_v3_fastqs.tar
  python3 ../split_by_barcode.py <lane 1 R1> <lane 1 R2> ../fastq 1
  python3 ../split_by_barcode.py <lane 2 R1> <lane 2 R2> ../fastq 2
  (needs pigz; e.g. in pipeline-bulkrnaseq-sm_env)

The previous inputs (shallow GBC353 samples, SLX-28283) are kept as inputs/*_SLX-28283.csv.

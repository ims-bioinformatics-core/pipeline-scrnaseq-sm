#!/bin/bash
# The account is not stored here: set it once in ~/.bashrc with  export SBATCH_ACCOUNT=<your account>
# (or pass it at submission: sbatch -A <your account> submit_slurm.sh)

#SBATCH -J pipeline-scrnaseq-sm         # Job name
#SBATCH --mem=78G                         # Max memory (icelake max 3370 MB/CPU -> 24 CPUs = ~79 GB; more mem = more CPUs billed)
#SBATCH --time=12:00:00                   # Max runtime (SL3 limit); resubmit to resume
#SBATCH --mail-type=ALL                   # Email notifications
#SBATCH -p icelake                        # Partition (76 cores, ~256 GB per node)
#SBATCH --cpus-per-task=24                # CPUs per task (FastQC: 12 samples x 2 threads; Cell Ranger: 1 sample x 24 cores / 64 GB)

# Usage: sbatch submit_slurm.sh   (from the pipeline directory)
bash run.sh

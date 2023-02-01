#!/bin/sh
#SBATCH --job-name=pb_distribution_test   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem-per-cpu=16G           # Memory per processor
#SBATCH --time=08:00:00             # Time limit hrs:min:sec
#SBATCH --output=logs/BS_test_%A-%a.out    # Standard output and error log
#SBATCH --array=2-41         # Array range

line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" sample_sheets/poisson_binomial_bootstrap_sample_sheet.txt)
#line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" sample_sheets/poisson_binomial_bootstrap_sample_sheet_2.txt)
#line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" poisson_binomial_bootstrap_sample_sheet_3.txt)

module load R/4.0.0

REP=$(echo $line | cut -f1 -d " ")
ALPHA=$(echo $line | cut -f2 -d " ")
BETA=$(echo $line | cut -f3 -d " ")
N=$(echo $line | cut -f4 -d " ")

Rscript pb-test-benchmark.R $REP $ALPHA $BETA $N


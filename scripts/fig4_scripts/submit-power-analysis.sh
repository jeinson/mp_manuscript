#!/bin/sh
#SBATCH --job-name=pb_distribution_test   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem-per-cpu=2G           # Memory per processor
#SBATCH --time=00:30:00             # Time limit hrs:min:sec
#SBATCH --output=logs/power_test_%A-%a.out    # Standard output and error log
#SBATCH --array=1-105          # Array range

line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" sample_sheets/pb_power_analysis_sample_sheet_fine.tsv)

module load R/4.0.0

ALPHA=$(echo $line | cut -f1 -d " ")
BETA=$(echo $line | cut -f2 -d " ")
SD_FROM_MEAN=$(echo $line | cut -f3 -d " ")

Rscript pb-power-test.R $ALPHA $BETA $SD_FROM_MEAN


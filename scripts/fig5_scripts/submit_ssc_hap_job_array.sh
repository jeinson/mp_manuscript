#!/bin/bash
#SBATCH --job-name=SSC_hap_config
#SBATCH --output=haplotype_configurations_v8_anno/ssc_hap_%A_%a.tsv
#SBATCH --error=error_logs/ssc_hap_stderrr_%A_%a.txt
#SBATCH --array=1-2380
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --mem=500M

# Print the task id.
#sample=$(head -${SLURM_ARRAY_TASK_ID} SSC_family_vcfs.list | tail -1) 
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p SSC_family_vcfs.list)

# Run the thing!
scripts/ssc_hap_configuration_executor.sh $sample

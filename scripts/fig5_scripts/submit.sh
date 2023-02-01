#!/bin/bash
#SBATCH --job-name=topmed_hap_config
#SBATCH --output=csnp_sqtl_haplotypes/topmed_hap_v8_anno_%A_%a.tsv
#SBATCH --error=error_logs/ssc_hap_stderrr_%A_%a.txt
#SBATCH --array=1-22
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G


GTEX_CODING_VCF="/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.coding.vcf.gz"
GTEX_ALL_VCF="/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.all.vcf.gz"
GTEX_EURO_VCF="/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_EUR_Analysis_Freeze.SHAPEIT2_phased.all.vcf.gz"
CADD="/gpfs/commons/groups/lappalainen_lab/data/cadd/v1.5/whole_genome_SNVs_inclAnno.tsv.gz"
ANC_ALLELE="/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/analysis/release_91_homo_sapiens.txt.gz"
basedir="/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/"
GTEX_UA_VCF="/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz"

# Select the VCF with the correct chromosome for this run
TOPMED_VCF="/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/TOPMED/topmed_vcfs/freeze.8.chr${SLURM_ARRAY_TASK_ID}.pass_only.phased_selected.vcf.gz"

source ~/python3_env/bin/activate
module load htslib/1.7

# Make sure to use the correct chromosome
GNOMAD="/gpfs/commons/datasets/gnomAD/3.0/gnomad.genomes.r3.0.sites.chr${SLURM_ARRAY_TASK_ID}.vcf.bgz"

scripts/haplotype_configuration_topmed.py \
--vcf_coding $TOPMED_VCF \
--vcf_all $TOPMED_VCF \
--cadd $CADD \
--gnomad $GNOMAD \
--anc_allele $ANC_ALLELE \
--genes /gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/multiple_exon_analysis/cross_tissue_top_sQTLs_and_secondary_exons_with_matching_signs_exons_labeled_by_coord.tsv \
--chrom "chr${SLURM_ARRAY_TASK_ID}"

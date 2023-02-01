#!/usr/bin/bash
# Setup some environment variables
GTEX_ALL_VCF="/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/ase/GTEx_Analysis_v8_phASER/phASER_GTEx_v8_merged.vcf.gz"
GTEX_EURO_VCF="/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_EUR_Analysis_Freeze.SHAPEIT2_phased.all.vcf.gz"
CADD="/gpfs/commons/groups/lappalainen_lab/data/cadd/v1.5/whole_genome_SNVs_inclAnno.tsv.gz"
ANC_ALLELE="/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/analysis/release_91_homo_sapiens.txt.gz"
basedir="/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/gtex_haplotypes_phaser"

# This was an argument before, but it makes things easier to hard code it. 
SQTLS="/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/multiple_exon_analysis/cross_tissue_top_sQTLs_and_secondary_exons_with_matching_signs_exons_labeled_by_coord.tsv"

set -eo pipefail

source /gpfs/commons/home/jeinson/python3_env/bin/activate
module load htslib/1.7


# check to make sure input is there
# if [[ ! -f ${basedir}/sQTL/combined_qtltools_results/${TISS}_combined_sQTLs.tsv ]]
# then
#     echo "Cannot find ${basedir}/sQTL/combined_qtltools_results/${TISS}_combined_sQTLs.tsv: Please check the path and try again"
#     exit
# fi

# if [[ ! -f ${basedir}/sQTL/indvs_per_tiss/${TISS}_indvs.txt ]]
# then 
#     echo "Cannot find ${basedir}/sQTL/indvs_per_tiss/${TISS}_indvs.txt: Please check the path and try again"
#     exit
# fi

a=$(seq 1 22)
for i in ${a[@]}; do
echo $i;

# Get the correct gnomad file
GNOMAD="/gpfs/commons/datasets/gnomAD/3.0/gnomad.genomes.r3.0.sites.chr$i.vcf.bgz"

sbatch --job-name hap_conf_${TISS}_$i -t 0:10:00 --mem=1G -e slurm_logs/cross_tiss_run_chr${i}.err \
	--output ${basedir}/haplotype_combos/secondary_exon_hap_combos_phaser_3_24_22_chr$i.tsv \
${basedir}/scripts/haplotype_configuration_gtex_phased_JE_multiple_exons.py \
--vcf_coding $GTEX_ALL_VCF \
--vcf_all $GTEX_ALL_VCF \
--cadd $CADD \
--gnomad $GNOMAD \
--anc_allele $ANC_ALLELE \
--genes $SQTLS \
--chrom $i

done

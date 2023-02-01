#!/usr/bin/bash

GTEX_CODING_VCF="/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.coding.vcf.gz"
GTEX_ALL_VCF="/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.all.vcf.gz"
GTEX_EURO_VCF="/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_EUR_Analysis_Freeze.SHAPEIT2_phased.all.vcf.gz"
CADD="/gpfs/commons/groups/lappalainen_lab/data/cadd/v1.5/whole_genome_SNVs_inclAnno.tsv.gz"
ANC_ALLELE="/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/analysis/release_91_homo_sapiens.txt.gz"
basedir="/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/ssc"
GTEX_UA_VCF="/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz"

# Get input file from the command line
SSC_SAMPLE=$1

source ~/python3_env/bin/activate
module load htslib/1.7

for i in $(seq 1 22); do
	GNOMAD="/gpfs/commons/datasets/gnomAD/3.0/gnomad.genomes.r3.0.sites.chr${i}.vcf.bgz"

	${basedir}/scripts/haplotype_configuration_gtex_multiple_exon_ssc.py \
	--vcf_coding $SSC_SAMPLE \
	--vcf_all $SSC_SAMPLE \
	--cadd $CADD \
	--gnomad $GNOMAD \
	--anc_allele $ANC_ALLELE \
	--genes ${basedir}/../sQTL_v8_anno/multiple_exon_analysis/cross_tissue_top_sQTLs_and_secondary_exons_with_matching_signs_exons_labeled_by_coord.tsv \
	--chrom $i 2> /dev/null;
done

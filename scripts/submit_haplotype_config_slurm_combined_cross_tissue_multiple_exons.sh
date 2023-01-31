#!/usr/bin/bash
# Setup some environment variables
GTEX_CODING_VCF="/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.coding.vcf.gz"
GTEX_ALL_VCF="/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.all.vcf.gz"
GTEX_EURO_VCF="/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_EUR_Analysis_Freeze.SHAPEIT2_phased.all.vcf.gz"
CADD="/gpfs/commons/groups/lappalainen_lab/data/cadd/v1.5/whole_genome_SNVs_inclAnno.tsv.gz"
ANC_ALLELE="/gpfs/commons/groups/lappalainen_lab/dglinos/projects/epistasis/analysis/release_91_homo_sapiens.txt.gz"
basedir="/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance"

set -eo pipefail

source /gpfs/commons/home/jeinson/python3_env/bin/activate
module load htslib/1.7

function Usage() {
    echo -e "\
Purpose: This script takes the path of combined sQTLs across tissues and runs haplotype calling
Usage: $(basename $0) -t path
Where: -t|--tiss is the relative path to the combined sQTL file
" >&2
    exit 1
}

# Janky argument parsing
[[ "$#" -lt 1 ]] && Usage
while [[ "$#" -gt 1 ]]; do
    case "$1" in
        -t|--tiss)
            TISS="$2"
            shift
            ;;
        *)
            Usage
            ;;
    esac
    shift
done

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
	--output ${basedir}/haplotype_combos/multiple_exons/secondary_exon_hap_combos_3_24_22_chr$i.tsv \
${basedir}/pipeline/scripts/haplotype_configuration_gtex_JE_multiple_exons.py \
--vcf_coding $GTEX_ALL_VCF \
--vcf_all $GTEX_ALL_VCF \
--cadd $CADD \
--gnomad $GNOMAD \
--anc_allele $ANC_ALLELE \
--genes $TISS \
--chrom $i

done

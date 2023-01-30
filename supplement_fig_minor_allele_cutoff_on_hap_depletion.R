#!/nfs/sw/R/R-4.0.0/bin/Rscript --vanilla

# This script performs haplotype enrichment analysis, but with increasingly 
# strict rare variant cutoffs. This visualizes how rare variants need to be
# to get a significant result. 



## This is the version of this script that lives in the mp_manuscript directory!!
rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/topmed/")
source("~/myPackages.R")
library(tompen)
library(ggplot2)
library(cowplot)
source("../mp_manuscript/tompen_utility_functions_manuscript.R")

# Read in all data
topmed_hap <- read_rds("topmed_haplotypes_v8_anno_full_with_details.rds")

#### Extract used cohorts from filtered set ####
# 6/15/22 Update: Filter to individuals in highly represented cohorts. 
# There were a few individuals from cohorts with mostly African Americans. 
# Getting approval for the whole cohort is a pain in the ass. 
used_indvs <- unique(topmed_hap$indv)
SRA_list <- read_tsv("cohort_info/selected_freeze8_sample_annot_2019-10-08.txt")
SRA_list_used_samples <- SRA_list %>% 
  filter(sample.id %in% used_indvs) %>%
  group_by(topmed_phs) %>% 
  summarise(n_samples = n())
topmed_studies <- readxl::read_xlsx("cohort_info/topmed_studies.xlsx")
SRA_list_used_samples$topmed_phs %in% 
  unique(topmed_studies$`TOPMed Accession # `)

# Perform the filter
n_samples_per_cohort <- deframe(SRA_list_used_samples)
cumulative_samples_per_cohort <- cumsum(sort(prop.table(n_samples_per_cohort), decreasing = T))

usable_cohorts <- names(cumulative_samples_per_cohort)[cumulative_samples_per_cohort < .95]
topmed_hap <- 
  topmed_hap %>%
  left_join(
    SRA_list %>% select(sample.id, topmed_phs),
    by = c("indv" = "sample.id")
  ) %>%
  filter(topmed_phs %in% usable_cohorts)

#### Now run the tests ####
f2 <- topmed_hap %>% filter(csnp_af_gnomad < 1e-2 & csnp_af_gnomad > 1e-3 & csnp_count > 1)
f3 <- topmed_hap %>% filter(csnp_af_gnomad < 1e-3 & csnp_af_gnomad > 1e-4 & csnp_count > 1)
f4 <- topmed_hap %>% filter(csnp_af_gnomad < 1e-4 & csnp_af_gnomad > 1e-5 & csnp_count > 1)
f5 <- topmed_hap %>% filter(csnp_af_gnomad < 1e-5 & csnp_count > 1)

# run the tests (this may take a while)
r2 <- run_enrichment(f2)
write_rds(r2, "cutoff_results/r2.rds")
r3 <- run_enrichment(f3)
write_rds(r3, "cutoff_results/r3.rds")
r4 <- run_enrichment(f4)
write_rds(r4, "cutoff_results/r4.rds")
r5 <- run_enrichment(f5)
write_rds(r5, "cutoff_results/r5.rds")

r2_comp <- run_comparison_test(f2)['bootstrap_p']
write_rds(r2_comp, "cutoff_results/r2_comp.rds")
r3_comp <- run_comparison_test(f3)['bootstrap_p']
write_rds(r3_comp, "cutoff_results/r3_comp.rds")
r4_comp <- run_comparison_test(f4)['bootstrap_p']
write_rds(r4_comp, "cutoff_results/r4_comp.rds")
r5_comp <- run_comparison_test(f5)['bootstrap_p']
write_rds(r5_comp, "cutoff_results/r5_comp.rds")

# Run 
# sbatch --mem 16g --time=10:00:00 supplement_fig_minor_allele_cutoff_on_hap_depletion.R

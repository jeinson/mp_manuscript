# Figure 7, which plots enrichment from SSC haplotypes. The pattern is the 
# same here.
# 
# We use the same types of filters used on GTEx and TOPMEd, for consistency
# across the paper at least. 

rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/ssc")
source("~/myPackages.R")
library("poibin")
source("../mp_manuscript/tompen_utility_functions_manuscript.R")
library(tompen)
library(ggplot2)

# Load the most recent haplotype configurations ----
# 
# Note that these have not been filtered yet!!! 
ssc_haps <- readRDS("haplotype_configurations_rds/SSC_sqtl_haplotypes_8_16_21.rds")

# Add information from sQTLs
top_sQTLs <- read_tsv("../sQTL_v8_anno/cross_tissue_top_sQTLs/top_sQTLs_MAF05.tsv")

ssc_haps <- 
  ssc_haps %>%
  left_join(
    top_sQTLs %>% select(top_pid, delta_psi, total_pid, mean_01_psi), 
    by = c("top_exon_coord" = "top_pid")
  )

# Add the expectation of beta, using gnomad allele frequencies. 
ssc_haps$exp_beta_empir <- 
  ifelse(ssc_haps$hom, 
         ssc_haps$esnp_gnomad_af^2 / 
           ((1-ssc_haps$esnp_gnomad_af)^2+ssc_haps$esnp_gnomad_af^2), 
         .5
  )

ssc_haps$exp_beta <- ssc_haps$exp_beta_empir

# Add a column for family
ssc_haps$family <- ssc_haps$indv %>% str_split("\\.") %>% map_chr(1)

# Split dataset by probands or unrelated parents
ssc_haps$haplotype <- factor(ssc_haps$haplotype, levels = names(beta_config_sqtl))
ssc_parents <- ssc_haps %>% filter(str_detect(indv, "mo|fa") )
ssc_proband <- ssc_haps %>% filter(str_detect(indv, "p1|p2"))
ssc_sibling <- ssc_haps %>% filter(str_detect(indv, "s1|s2"))

### Draw summary plots before filtering. 
n_haps_per_gene <- table(ssc_parents$gene)

save_plot("fig7S1_N_haps_per_gene_ssc_prefilter.svg", width = 5, height = 3)
boxplot(as.numeric(n_haps_per_gene), horizontal = T, log = 'x', 
        xlab = "Number of haplotypes per gene in SSC parents")
dev.off()

# How to properly filter these haplotype calls to focus on actual rare variants
parent_csnp_table <- table(ssc_parents$csnp)
rarest_variants <- names(parent_csnp_table[which(parent_csnp_table <= 3)])

# Filter out genes that have an UNUSUALLY high number of haplotypes across
# parents. 
outlier_genes <- flag_outliers(table(ssc_parents$gene))
usable_ssc_genes <- true_names(!outlier_genes)

ssc_parents <- filter(ssc_parents, csnp %in% rarest_variants & gene %in% usable_ssc_genes)
ssc_proband <- filter(ssc_proband, csnp %in% rarest_variants & gene %in% usable_ssc_genes)
ssc_sibling <- filter(ssc_sibling, csnp %in% rarest_variants & gene %in% usable_ssc_genes)

save_plot("fig7S2_N_haps_per_indv_ssc_parents.svg", width = 4, height = 4)
hist(as.numeric(table(ssc_parents$indv)), col = "chartreuse3", breaks = 20)
box()
dev.off()

# First check depletion across all variants and genes
poison_binomial_test(ssc_parents)

# Now run the enrichment and plot the plots!
library(tompen)
ssc_parents_result <- run_enrichment(ssc_parents)
ssc_proband_result <- run_enrichment(ssc_proband)
ssc_sibling_result <- run_enrichment(ssc_sibling)

ssc_parents_comp <- run_comparison_test(ssc_parents)['bootstrap_p']
ssc_proband_comp <- run_comparison_test(ssc_proband)['bootstrap_p']
ssc_sibling_comp <- run_comparison_test(ssc_sibling)['bootstrap_p']

library(cowplot)
save_plot("fig7A_ssc_overall_enrichment.svg", width = 5, height = 5)
plot_grid(
  plot_epsilon_plot(ssc_parents_result, ssc_parents_comp, 
                    pval_pos = .09, n_pos = .12) + 
    ggtitle("Haplotypes depletion in SSC Parents"),
  plot_epsilon_plot(ssc_proband_result, ssc_proband_comp,
                    pval_pos = .09, n_pos = .12) + 
    ggtitle("Haplotype depletion in SSC Probands"),
  plot_epsilon_plot(ssc_sibling_result, ssc_sibling_comp, 
                    pval_pos = .09, n_pos = .12) + 
    ggtitle("Haplotype depletion in SSC Siblings"),
  nrow = 3
)
dev.off()

# Check variants shared between probands and siblings, like Stephan did in 
# his paper. 
proband_csnps <- unique(ssc_proband$csnp)
sibling_csnps <- unique(ssc_sibling$csnp)

proband_unique_csnps <- setdiff(proband_csnps, sibling_csnps)
sibling_unique_csnps <- setdiff(sibling_csnps, proband_csnps)
shared_csnps <- intersect(proband_csnps, sibling_csnps)

run_enrichment(ssc_proband %>% filter(csnp %in% proband_unique_csnps))
run_enrichment(ssc_sibling %>% filter(csnp %in% sibling_unique_csnps))

run_enrichment(ssc_proband %>% filter(csnp %in% shared_csnps))
run_enrichment(ssc_sibling %>% filter(csnp %in% shared_csnps))

# Compare depletion of these variants to topmed...?
topmed_hap <- read_rds("../topmed/topmed_haplotypes_v8_anno_full_with_details.rds")
tm_to_ssd_snp_convert <- function(x){
  x = str_split(x, pattern = "_")
  map_chr(x, ~ paste0(.x[1], ":", .x[2]))
}

topmed_hap <- 
  topmed_hap %>%
  filter(csnp_count <= 10 & csnp_count > 1) %>%
  mutate(ssc_snp_id = tm_to_ssd_snp_convert(csnp)) %>%
  mutate(exp_beta = exp_beta_empir)

ssc_proband_pb_unique <- ssc_proband %>% filter(csnp %in% proband_unique_csnps)
topmed_pb_unique <- topmed_hap %>% filter(ssc_snp_id %in% proband_unique_csnps)

tm_r0 <- run_enrichment(topmed_pb_unique)
ssc_r0 <- run_enrichment(ssc_proband_pb_unique)
tm_r0_comp <- run_comparison_test(topmed_pb_unique)
ssc_r0_comp <- run_comparison_test(ssc_proband_pb_unique)

save_plot("fig7B_topmed_SSC_comparison_all_genes.svg", width = 5, height = 3)
plot_grid(
  plot_epsilon_plot(x = ssc_r0, comp_pval = ssc_r0_comp['bootstrap_p'], 
                    pval_pos = .04, n_pos = .055) + 
    ggtitle("SSC proband haplotypes with cSNPs only in SSC probands"),
  plot_epsilon_plot(x = tm_r0, comp_pval = tm_r0_comp['bootstrap_p'], 
                    pval_pos = .04, n_pos = .055) + 
    ggtitle("TOPMEd haplotypes with cSNPs only in SSC probands"),
  nrow = 2
)
dev.off()

#### Check SFARI Autism Genes ####
# sfari_genes <- read_csv("SFARI-Gene_genes_08-07-2020release_08-14-2020export.csv")
sfari_genes <- read_csv("SFARI-Gene_genes_09-02-2021release_10-25-2021export.csv")
sfari_genes <- sfari_genes %>% filter(!is.na(`ensembl-id`))
sfari_genes$`gene-score`[is.na(sfari_genes$`gene-score`)] <- 4
l3_genes <- with(sfari_genes, `ensembl-id`[`gene-score` == 3])
l2_genes <- with(sfari_genes, `ensembl-id`[`gene-score` == 2])
l1_genes <- with(sfari_genes, `ensembl-id`[`gene-score` == 1])

# How many of each category do we actually have haps for?
sum(l3_genes %in% ssc_haps$gene) / length(l3_genes)
sum(l2_genes %in% ssc_haps$gene) / length(l2_genes)
sum(l1_genes %in% ssc_haps$gene) / length(l1_genes)

proband_sfari_gene_result <- ssc_proband %>%
  filter(gene %in% sfari_genes$`ensembl-id`) %>%
  run_enrichment
sibling_sfari_gene_result <- ssc_sibling %>%
  filter(gene %in% sfari_genes$`ensembl-id`) %>%
  run_enrichment
parents_sfari_gene_result <- ssc_parents %>%
  filter(gene %in% sfari_genes$`ensembl-id`) %>%
  run_enrichment

proband_sfari_gene_comp_result <- ssc_proband %>%
  filter(gene %in% sfari_genes$`ensembl-id`) %>%
  run_comparison_test() %>% .["bootstrap_p"]
sibling_sfari_gene_comp_result <- ssc_sibling %>%
  filter(gene %in% sfari_genes$`ensembl-id`) %>%
  run_comparison_test() %>% .["bootstrap_p"]
parents_sfari_gene_comp_result <- ssc_parents %>%
  filter(gene %in% sfari_genes$`ensembl-id`) %>%
  run_comparison_test() %>% .["bootstrap_p"]


# Plot enrichment for just SSC genes
plot_grid(
  plot_epsilon_plot(parents_sfari_gene_result, parents_sfari_gene_comp_result,
                    pval_pos = .09, n_pos = .12) + 
    ggtitle("SSC Parents, ASD relevant genes"),
  plot_epsilon_plot(proband_sfari_gene_result, proband_sfari_gene_comp_result, 
                    pval_pos = .09, n_pos = .12) +
    ggtitle("SSC Probands, ASD relevant genes"),
  plot_epsilon_plot(sibling_sfari_gene_result, sibling_sfari_gene_comp_result,
                    pval_pos = .09, n_pos = .12) + 
    ggtitle("SSC Siblings, ASD relevant genes"),
  
  nrow = 3
)

# Check these but with variants unique to probands. 
proband_sfari_gene_result_pb_only <- ssc_proband %>%
  filter(gene %in% sfari_genes$`ensembl-id` & csnp %in% proband_unique_csnps) %>%
  run_enrichment
proband_sfari_gene_result_pb_only_comp <- ssc_proband %>%
  filter(gene %in% sfari_genes$`ensembl-id` & csnp %in% proband_unique_csnps) %>%
  run_comparison_test() %>% .["bootstrap_p"]

topmed_sfari_gene_result_pb_only <- topmed_hap %>%
  filter(gene %in% sfari_genes$`ensembl-id` & ssc_snp_id %in% proband_unique_csnps) %>%
  run_enrichment
topmed_sfari_gene_result_pb_only_comp <- topmed_hap %>%
  filter(gene %in% sfari_genes$`ensembl-id` & ssc_snp_id %in% proband_unique_csnps) %>%
  run_comparison_test() %>% .["bootstrap_p"]

save_plot("fig7C_topmed_SSC_comparison_ASD_genes.svg", width = 5, height = 3)
plot_grid(
  plot_epsilon_plot(proband_sfari_gene_result_pb_only, proband_sfari_gene_result_pb_only_comp, 
                    pval_pos = .12, n_pos = .14) + 
    ggtitle("SSC probands, ASD genes, cSNPs only in SSC probands"), 
  plot_epsilon_plot(topmed_sfari_gene_result_pb_only, topmed_sfari_gene_result_pb_only_comp, 
                    pval_pos = .12, n_pos = .14) + 
    ggtitle("TOPMED donors, ASD genes, cSNPs only in SSC probands"),
  nrow = 2
)
dev.off()


#### Redo these analyses, but in a more organized way for the purpose of the paper
# 1) All genes, SSC probands vs. tomed (del + non-del)
ssc_proband <- filter(ssc_proband, csnp %in% rarest_variants & gene %in% usable_ssc_genes)
ssc_proband_result <- run_enrichment(ssc_proband)
ssc_proband_comp <- run_comparison_test(ssc_proband)['bootstrap_p']

topmed_ssc_vars_genes <- filter(topmed_hap, ssc_snp_id %in% rarest_variants & gene %in% usable_ssc_genes)

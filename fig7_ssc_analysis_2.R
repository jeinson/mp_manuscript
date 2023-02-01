# Produce figures related to SSC data used in the main analysis. 

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
# 
# Load Data and apply filters ----
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

# How to properly filter these haplotype calls to focus on actual rare variants
parent_csnp_table <- table(ssc_parents$csnp)
rarest_variants <- names(parent_csnp_table[which(parent_csnp_table <= 3)])

# Filter out genes that have an UNUSUALLY high number of haplotypes across
# parents. 
outlier_genes <- flag_outliers(table(ssc_parents$gene))
usable_ssc_genes <- true_names(!outlier_genes)

# Apply the filters
ssc_parents <- filter(ssc_parents, csnp %in% rarest_variants & gene %in% usable_ssc_genes)
ssc_proband <- filter(ssc_proband, csnp %in% rarest_variants & gene %in% usable_ssc_genes)
ssc_sibling <- filter(ssc_sibling, csnp %in% rarest_variants & gene %in% usable_ssc_genes)

#### Make the 5 summary plots for the supplemenmt ####
## High inclusion allele frequency 
hi_af <- 
  ssc_parents %>%
  select(gene, esnp_gnomad_af) %>%
  distinct %>%
  ggplot(aes(esnp_gnomad_af)) +
  geom_histogram(fill = "lightgrey", color = "black") +
  xlab("Frequency of higher exon inclusion allele (gnomad)") +
  ylab("ψQTLs") + 
  gtex_v8_figure_theme() 

## Haplotypes per individual
x <- ssc_parents %>%
  group_by(indv) %>%
  summarize(n_haps = n()) %>%
  ungroup 

hap_per_indv <- 
  ggplot(x, aes(n_haps)) + 
  geom_histogram(fill = "lightgrey", color = "black", bins = 20) + 
  xlim(0, 50) +
  xlab("N Haplotypes per donor") + 
  ylab("Donors") +
  gtex_v8_figure_theme()

## Haplotypes per gene
x <- ssc_parents %>%
  mutate(del = ifelse(CADD > 15, "Pathogenic", "Benign")) %>%
  group_by(gene, del) %>%
  summarize(n_haps = n()) %>%
  ungroup 

hap_p_gene <- 
  ggplot(x, aes(log10(n_haps), fill = del)) + 
  geom_histogram(data = subset(x, del == "Pathogenic"), fill = "red", color = "black", alpha = .2, binwidth = .25) + 
  geom_histogram(data = subset(x, del == "Benign"), fill = "blue", color = "black", alpha = .2, binwidth = .25) + 
  xlab("log10(Haplotypes per Gene)") + 
  ylab("Genes") +
  xlim(0, 3.5) +
  gtex_v8_figure_theme()

## Coding variant allele frequencies. 
x <- ssc_parents %>%
  mutate(del = ifelse(CADD > 15, "Pathogenic", "Benign")) %>%
  select(csnp, csnp_af_gnomad, del) %>%
  distinct %>%
  ungroup 

cv_af <- ggplot(x, aes(csnp_af_gnomad, fill = del)) + 
  geom_histogram(data = subset(x, del == "Pathogenic"), fill = "red", color = "black", alpha = .2) + 
  geom_histogram(data = subset(x, del == "Benign"), fill = "blue", color = "black", alpha = .2) + 
  xlim(0, .0005) +
  xlab("Coding variant frequency in gnomad") + 
  ylab("Genes") +
  gtex_v8_figure_theme()

## Plot everything together
library(cowplot)
save_plot("figS3_ssc_parent_haplotype_summaries.svg", width = 10, height = 2)
plot_grid(hi_af, hap_per_indv, hap_p_gene, cv_af, nrow = 1)
dev.off()



# Load TOPMEd data ----
topmed_hap <- read_rds("../topmed/topmed_haplotypes_v8_anno_full_with_details.rds")
tm_to_ssd_snp_convert <- function(x){
  x = str_split(x, pattern = "_")
  map_chr(x, ~ paste0(.x[1], ":", .x[2]))
}

SRA_list <- read_tsv(here("topmed/cohort_info/selected_freeze8_sample_annot_2019-10-08.txt"))
used_cohorts <- read_lines(here("topmed/cohort_info/19_topmed_cohorts_in_Einson_2022.txt"))

topmed_hap %<>%
  left_join(
    SRA_list %>% select(sample.id, topmed_phs), 
    by = c("indv" = "sample.id")
  ) %>% 
  filter(topmed_phs %in% used_cohorts)

topmed_hap <- 
  topmed_hap %>%
  filter(csnp_count <= 10 & csnp_count > 1) %>%
  mutate(ssc_snp_id = tm_to_ssd_snp_convert(csnp)) %>%
  mutate(exp_beta = exp_beta_empir)

# Apply the same filters used on the SSC data
topmed_hap <- filter(topmed_hap, ssc_snp_id %in% rarest_variants & gene %in% usable_ssc_genes)

# #### 1) All genes in SSC vs. TOPMed ####
# 
# # First test for depletion in SSC Parents
# p1_parent <- run_enrichment(ssc_parents)
# p1_parent_comp <- run_comparison_test(ssc_parents)
# 
# # Then do probands, and TOPMEd limited to the same set
# p1_sp <- run_enrichment(ssc_proband)
# p1_sp_comp <- run_comparison_test(ssc_proband)  
# p1_tm <- run_enrichment(topmed_hap)
# p1_tm_comp <- run_comparison_test(topmed_hap)
# 
# p1_comp <- 
#   bootstrap_comparison_test(
#     ssc_proband %>% filter(CADD > 15), 
#     topmed_hap %>% filter(CADD > 15)
#   )
# 
# # Check in variants between the two
# prop.table(table(unique(topmed_hap$ssc_snp_id) %in% unique(ssc_proband$csnp)))
# 
# #### 2) Variants unique to SSC Probands (compared to siblings) ####
# proband_csnps <- unique(ssc_proband$csnp)
# sibling_csnps <- unique(ssc_sibling$csnp)
# 
# proband_unique_csnps <- setdiff(proband_csnps, sibling_csnps)
# sibling_unique_csnps <- setdiff(sibling_csnps, proband_csnps)
# shared_csnps <- intersect(proband_csnps, sibling_csnps)
# 
# p2_sp <- ssc_proband %>% filter(csnp %in% proband_unique_csnps) %>% run_enrichment()
# p2_sp_comp <- ssc_proband %>% filter(csnp %in% proband_unique_csnps) %>% run_comparison_test()
# p2_tm <- topmed_hap %>% filter(ssc_snp_id %in% proband_unique_csnps) %>% run_enrichment()
# p2_tm_comp <- topmed_hap %>% filter(ssc_snp_id %in% proband_unique_csnps) %>% run_comparison_test()
# 
# p2_comp <- bootstrap_comparison_test(
#   ssc_proband %>% filter(csnp %in% proband_unique_csnps & CADD > 15), 
#   topmed_hap %>% filter(ssc_snp_id %in% proband_unique_csnps & CADD > 15)
# )
# 
# # Check in variants between the two
# prop.table(table(
#   unique(filter(ssc_proband, csnp %in% proband_unique_csnps)$csnp) %in% 
#     unique(filter(topmed_hap, ssc_snp_id %in% proband_unique_csnps)$ssc_snp_id)))
# 
# #### 3) Variants shared between SSC Probands and Unaffected Siblings ####
# p3_sp <- ssc_proband %>% filter(csnp %in% shared_csnps) %>% run_enrichment()
# p3_sp_comp <- ssc_proband %>% filter(csnp %in% shared_csnps) %>% run_comparison_test()
# p3_tm <- topmed_hap %>% filter(ssc_snp_id %in% shared_csnps) %>% run_enrichment()
# p3_tm_comp <- topmed_hap %>% filter(ssc_snp_id %in% shared_csnps) %>% run_comparison_test()
# 
# p3_comp <- bootstrap_comparison_test(
#   ssc_proband %>% filter(csnp %in% shared_csnps & CADD > 15),
#   topmed_hap %>% filter(ssc_snp_id %in% shared_csnps & CADD > 15)
# )
# 
# # Check in variants between the two
# prop.table(table(
#   unique(filter(ssc_proband, csnp %in% shared_csnps)$csnp) %in% 
#     unique(filter(topmed_hap, ssc_snp_id %in% shared_csnps)$ssc_snp_id)))

#### 4) Variants in ASD specific genes ----
sfari_genes <- read_csv("SFARI-Gene_genes_09-02-2021release_10-25-2021export.csv")
sfari_genes <- sfari_genes %>% filter(!is.na(`ensembl-id`))
# 
# p4_sp <- ssc_proband %>% 
#   filter(csnp %in% proband_csnps & gene %in% sfari_genes$`ensembl-id`) %>% 
#   run_enrichment()
# p4_sp_comp <- ssc_proband %>% 
#   filter(csnp %in% proband_csnps & gene %in% sfari_genes$`ensembl-id`) %>% 
#   run_comparison_test()
# p4_tm <- topmed_hap %>%
#   filter(ssc_snp_id %in% shared_csnps & gene %in% sfari_genes$`ensembl-id`) %>%
#   run_enrichment()
# p4_tm_comp <- topmed_hap %>%
#   filter(ssc_snp_id %in% shared_csnps & gene %in% sfari_genes$`ensembl-id`) %>%
#   run_comparison_test()
# 
# p4_comp <- bootstrap_comparison_test(
#   ssc_proband %>%
#     filter(csnp %in% proband_csnps & gene %in% sfari_genes$`ensembl-id` & CADD > 15),
#   topmed_hap %>%
#     filter(ssc_snp_id %in% shared_csnps & gene %in% sfari_genes$`ensembl-id` & CADD > 15)
# )

## Check if ASD genes are more constrained. 
gene_constraints <- read_tsv("/gpfs/commons/groups/lappalainen_lab/jeinson/data/gnomad.v2.1.1.lof_metrics.by_gene.txt")
gene_map <- read_tsv("/gpfs/commons/groups/lappalainen_lab/jmorris/ensembl/ensemblKey.txt")
gene_constraints %<>%
  rename("name" = "gene") %>%
  left_join(gene_map) %>%
  select(gene, name, oe_lof_upper)

sfari_gene_constraint <- left_join(gene_constraints, sfari_genes, 
                                   by = c("gene" = "ensembl-id"))
sfari_gene_constraint$is_asd_gene <- !is.na(sfari_gene_constraint$`gene-name`)

asd_gene_LOEUF_density <- 
  ggplot(sfari_gene_constraint, aes(oe_lof_upper, fill = is_asd_gene)) + 
  geom_density(alpha = .5) + 
  xlab("LOEUF") + ylab("density (genes)") +
  xlim(0, 3) + 
  gtex_v8_figure_theme() + 
  theme(legend.position = c(.75,.75))

t.test(oe_lof_upper ~ is_asd_gene, data = sfari_gene_constraint)

### Draw up the plots

# Draw up the enrichment in Parents, which will be used to compare the other
# two groups to. 
pval_pos = .023
n_pos = 0.028

save_plot("fig5A_ssc_parent_haplotype_enrichment.svg", width = 5, height = 1.25)
plot_epsilon_plot(p1_parent, pval_pos = pval_pos, n_pos = n_pos,comp_pval = p1_parent_comp['bootstrap_p'])
dev.off()

# Now do the other set of analyses
# pval_pos = .1
# n_pos = .13
# save_plot("fig7_SSC_haplotype_depletion_cross.svg", width = 5, height = 13.5)
# plot_grid(
#   plot_epsilon_plot(p1_parent, p1_parent_comp['bootstrap_p'], pval_pos, n_pos) + ggtitle("p1_parent"),
#   
#   plot_epsilon_plot(p1_sp, p1_sp_comp['bootstrap_p'], pval_pos, n_pos) + ggtitle("p1_sp"),
#   plot_epsilon_plot(p1_tm, p1_tm_comp['bootstrap_p'], pval_pos, n_pos) + ggtitle("p1_tm"),
#   
#   plot_epsilon_plot(p2_sp, p2_sp_comp['bootstrap_p'], pval_pos, n_pos) + ggtitle("p2_sp"),
#   plot_epsilon_plot(p2_tm, p2_tm_comp['bootstrap_p'], pval_pos, n_pos) + ggtitle("p2_tm"),
#   
#   plot_epsilon_plot(p3_sp, p3_sp_comp['bootstrap_p'], pval_pos, n_pos) + ggtitle("p3_sp"),
#   plot_epsilon_plot(p3_tm, p3_tm_comp['bootstrap_p'], pval_pos, n_pos) + ggtitle("p3_tm"),
#   
#   plot_epsilon_plot(p4_sp, p4_sp_comp['bootstrap_p'], pval_pos, n_pos) + ggtitle("p4_sp"),
#   plot_epsilon_plot(p4_tm, p4_tm_comp['bootstrap_p'], pval_pos, n_pos) + ggtitle("p4_tm"),
#   
#   ncol = 1
# )
# dev.off()
# 
# rbind(p1_comp, p2_comp, p3_comp, p4_comp)

save_plot("figS9B_constrained_gene_loeuf.svg", width = 3, height = 2)
asd_gene_LOEUF_density
dev.off()

# Now it's time to re-do this analysis, because apparently it got 
# deleted somewhere along the line. 

sfari_parent  <- ssc_parents %>%
  filter(gene %in% sfari_genes$`ensembl-id`) %>% 
  run_enrichment()
sfari_parent_comp <- ssc_parents %>% 
  filter(gene %in% sfari_genes$`ensembl-id`) %>% 
  run_comparison_test()

sfari_proband <- ssc_proband %>% 
  filter(gene %in% sfari_genes$`ensembl-id`) %>% 
  run_enrichment()
sfari_proband_comp <- ssc_proband %>% 
  filter(gene %in% sfari_genes$`ensembl-id`) %>% 
  run_comparison_test()

sfari_sibling <- ssc_sibling %>% 
  filter(gene %in% sfari_genes$`ensembl-id`) %>% 
  run_enrichment()
sfari_sibling_comp <- ssc_sibling %>% 
  filter(gene %in% sfari_genes$`ensembl-id`) %>% 
  run_comparison_test()

sfari_tm <- topmed_hap %>%
  filter(gene %in% sfari_genes$`ensembl-id`) %>%
  run_enrichment()
sfari_tm_comp <- topmed_hap %>%
  filter(gene %in% sfari_genes$`ensembl-id`) %>%
  run_comparison_test()

library(cowplot)
pval_pos = .12
n_pos = 0.1
fig7_topmed_SSC_comp_ASD_genes <- 
  plot_grid(
    plot_epsilon_plot(sfari_parent, sfari_parent_comp['bootstrap_p'], n_pos, pval_pos) + ggtitle("Parents"),
    plot_epsilon_plot(sfari_proband, sfari_proband_comp['bootstrap_p'], n_pos, pval_pos) + ggtitle("Probands"),
    plot_epsilon_plot(sfari_sibling, sfari_sibling_comp['bootstrap_p'], n_pos, pval_pos) + ggtitle("Siblings"),
    plot_epsilon_plot(sfari_tm, sfari_tm_comp['bootstrap_p'], n_pos, pval_pos) + ggtitle("TOPMed"),
    
    nrow = 4
)

save_plot("fig7_SSC_haplotype_depletion_cross.svg", width = 5, height = 6)
fig7_topmed_SSC_comp_ASD_genes
dev.off()

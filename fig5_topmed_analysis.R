# This script performs the same analysis as in full_topmed_analysis.R, with more
# emphasis on code sharing and figure tweaking. 

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

for_tuuli <- inner_join(topmed_studies, SRA_list_used_samples, 
                        by = c(`TOPMed Accession # ` = "topmed_phs")) %>%
  group_by(`TOPMed Accession # `) %>%
  slice(1) %>%
  ungroup
write_tsv(for_tuuli, "cohort_info/topmed_studies_used_in_Einson_2022.tsv")

# Perform the filter
n_samples_per_cohort <- deframe(SRA_list_used_samples)
cumulative_samples_per_cohort <- cumsum(sort(prop.table(n_samples_per_cohort), decreasing = T))
plot(cumulative_samples_per_cohort, main = "Cumulative percentage of European samples across TOPMed cohorts")
abline(h = .95, col = "red")
grid()

usable_cohorts <- names(cumulative_samples_per_cohort)[cumulative_samples_per_cohort < .95]
topmed_hap <- 
  topmed_hap %>%
  left_join(
    SRA_list %>% select(sample.id, topmed_phs),
    by = c("indv" = "sample.id")
  ) %>%
  filter(topmed_phs %in% usable_cohorts)

# Impose a stricter minor allele frequency filter. I only want variants that
# appear 10 or fewer times in the whole dataset. 
topmed_hap <- filter(topmed_hap, csnp_count <= 10 & csnp_count > 1)

# Final filtering procedure: Filter out genes with a large number of haplotypes 
# across individuals (not across eexons!)
haps_per_gene <- table(topmed_hap$gene)
usable_genes <- true_names(!flag_outliers(haps_per_gene))
topmed_hap <- topmed_hap %>% filter(gene %in% usable_genes)

#### Make the 5 summary plots for the supplemenmt ####
## High inclusion allele frequency 
hi_af <- 
  topmed_hap %>%
  select(gene, recom_esnp_af) %>%
  distinct %>%
  ggplot(aes(recom_esnp_af)) +
  geom_histogram(fill = "lightgrey", color = "black") +
  xlab("Frequency of higher exon inclusion allele") +
  ylab("ψQTLs") + 
  gtex_v8_figure_theme() 

## Haplotypes per individual
x <- topmed_hap %>%
  group_by(indv) %>%
  summarize(n_haps = n()) %>%
  ungroup 

hap_per_indv <- 
  ggplot(x, aes(n_haps)) + 
  geom_histogram(fill = "lightgrey", color = "black") + 
  xlim(0, 50) +
  xlab("N Haplotypes per donor") + 
  ylab("Donors") +
  gtex_v8_figure_theme()

## Haplotypes per gene
x <- topmed_hap %>%
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
x <- topmed_hap %>%
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
save_plot("figS3_topmed_haplotype_summaries.svg", width = 10, height = 2)
plot_grid(hi_af, hap_per_indv, hap_p_gene, cv_af, nrow = 1)
dev.off()

# Constraint metrics, in quantiles
loeuf_quantiles <- select(topmed_hap, gene, oe_lof) %>% distinct %$% 
  quantile(oe_lof, na.rm = T)

dpsi_quantiles <- 
  select(topmed_hap, gene, delta_psi) %>%
  distinct %$%
  quantile(delta_psi)

#### Summary Plots ####
topmed_hap <- 
  within(topmed_hap, 
         effect <-
           factor(effect,
                  levels = names(sort(table(effect), decreasing = T)))
  )

topmed_hap <- 
  topmed_hap %>%
  mutate(deleterious = ifelse(CADD > 15, "deleterious", "non-deleterious"))

save_plot("fig5S2_RV_annotation_counts.svg", width = 4, height = 3)
ggplot(topmed_hap, aes(x = effect, fill = deleterious)) + 
  geom_bar() + 
  scale_fill_brewer("CADD Designation", palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("Rare variant counts") +
  gtex_v8_figure_theme() + 
  theme(legend.position = c(.7, .7))
dev.off()

### MP by filters ----
# Add this to keep the thing from breaking
topmed_hap$exp_beta <- topmed_hap$exp_beta_empir

# Run depletion across EVERYTHING
mean(topmed_hap$beta - topmed_hap$exp_beta)
poibin::ppoibin(sum(topmed_hap$beta), topmed_hap$exp_beta_empir, method = "RNA") * 2

# Separate the haplotypes based on filters
f1 <- filter(topmed_hap, delta_psi > dpsi_quantiles[4])
f2 <- filter(topmed_hap, oe_lof < loeuf_quantiles[2])
f3 <- filter(topmed_hap, delta_psi > dpsi_quantiles[4] &
               oe_lof < loeuf_quantiles[2])

#### Run analysis and plot ----
# Actually run the thing! The first filters used
r0 <- run_enrichment(topmed_hap)
r1 <- run_enrichment(f1)
r2 <- run_enrichment(f2)
r3 <- run_enrichment(f3)

# Now run comparison analyses
r0_comp <- run_comparison_test(topmed_hap)['bootstrap_p']
r1_comp <- run_comparison_test(f1)['bootstrap_p']
r2_comp <- run_comparison_test(f2)['bootstrap_p']
r3_comp <- run_comparison_test(f3)['bootstrap_p']

# Plot the plots
pval_pos = .023
n_pos = 0.028

#### first round of filtering plot ####
save_plot("fig5A_topmed_haps_by_3_filters.svg", width = 5, height = 5)
plot_grid(
  plot_epsilon_plot(r0, comp_pval = r0_comp, pval_pos = pval_pos, n_pos = n_pos),
  plot_epsilon_plot(r1, comp_pval = r1_comp, pval_pos = pval_pos, n_pos = n_pos),
  plot_epsilon_plot(r2, comp_pval = r2_comp, pval_pos = pval_pos, n_pos = n_pos),
  plot_epsilon_plot(r3, comp_pval = r3_comp, pval_pos = pval_pos, n_pos = n_pos), 
  ncol = 1
)
dev.off()
#####

#### Homs vs. Hets. ####
r0_het <- run_enrichment(topmed_hap %>% filter(!hom))
r0_het_comp <- run_comparison_test(topmed_hap %>% filter(!hom))
r0_hom <- run_enrichment(topmed_hap %>% filter(hom))
r0_hom_comp <- run_comparison_test(topmed_hap %>% filter(hom))

save_plot("fig5S4_haps_by_hom_or_het.svg", width = 5, height = 3)
plot_grid(
  plot_epsilon_plot(x = r0_het, comp_pval = r0_het_comp["bootstrap_p"]) + ggtitle("Heterozygotes"),
  plot_epsilon_plot(x = r0_hom, comp_pval = r0_hom_comp["bootstrap_p"]) + ggtitle("Homozygotes"), 
  nrow = 2
)
dev.off()

### Derived allele effect direction ----
topmed_hap$der_allele_hl <- 
  with(topmed_hap,
       # This is based on the ancestral allele, 
       # so flip it since we're interested in the derived allele
       ifelse(esnp_anc_allele == esnp_hi_inc_allele,
              "lower",
              ifelse(esnp_anc_allele == esnp_low_inc_allele,
                     "higher", NA
              )
       )
  )

r0_lower <- run_enrichment(topmed_hap %>% filter(der_allele_hl == "lower"))
r0_lower_comp <- run_comparison_test(topmed_hap %>% filter(der_allele_hl == "lower"))['bootstrap_p']
r0_higher <- run_enrichment(topmed_hap %>% filter(der_allele_hl == "higher"))
r0_higher_comp <- run_comparison_test(topmed_hap %>% filter(der_allele_hl == "higher"))['bootstrap_p']

save_plot("fig5S5_derived_allele_direction.svg", width = 5, height = 3)
plot_grid(
  plot_epsilon_plot(x = r0_lower, comp_pval = r0_lower_comp) + ggtitle("sQTL derived allele decreases exon inclusion"),
  plot_epsilon_plot(x = r0_higher, comp_pval = r0_higher_comp) + ggtitle("sQTL derived allele increases exon inclusion"), 
  nrow = 2
)
dev.off()

### Overall exon inclusion enrichment ----
r0_mostly_spliced <- run_enrichment(topmed_hap %>% filter(mean_01_psi < .4))
r0_mostly_spliced_comparison <- 
  run_comparison_test(topmed_hap %>% filter(mean_01_psi < .4))['bootstrap_p']
r0_mostly_included <- run_enrichment(topmed_hap %>% filter(mean_01_psi > .6))
r0_mostly_included_comparison <- 
  run_comparison_test(topmed_hap %>% filter(mean_01_psi > .6))['bootstrap_p']

#### overall inclusion figure ####
save_plot("fig5S6_overall_inclusion.svg", width = 5, height = 3)
plot_grid(
  plot_epsilon_plot(r0_mostly_spliced, comp_pval = r0_mostly_spliced_comparison) + 
    ggtitle("sQTLs affecting overall lowly included exons"),
  plot_epsilon_plot(r0_mostly_included, comp_pval = r0_mostly_included_comparison) + 
    ggtitle("sQTL affecting overall highly included exons"), 
  nrow = 2
)
dev.off()
#####

# As a side note: 
save_plot("fig5S7_snp_cadd_density.svg", width = 5, height = 4)
topmed_hap %>% 
  select(csnp, mean_01_psi, CADD) %>% 
  distinct %>% 
  ggplot(aes(mean_01_psi, CADD)) + 
  #geom_point() +
  geom_density_2d(col = "black") +
  geom_density_2d_filled(alpha = .5) + 
  geom_smooth() + 
  gtex_v8_figure_theme()
dev.off()

### MP as a function of sQTL effect size ----
dpsi_quant <- 
  topmed_hap %>% select(gene, delta_psi) %>% distinct %$% quantile(delta_psi)
c1 = round(dpsi_quant[2], 3)
c2 = round(dpsi_quant[3], 3)
c3 = round(dpsi_quant[4], 3)
r0_q1_dpsi <- run_enrichment(topmed_hap %>% filter(delta_psi < c1))
r0_q2_dpsi <- run_enrichment(topmed_hap %>% filter(between(delta_psi, c1, c2)))
r0_q3_dpsi <- run_enrichment(topmed_hap %>% filter(between(delta_psi, c2, c3)))
r0_q4_dpsi <- run_enrichment(topmed_hap %>% filter(delta_psi > c3))
# comparisons
r0_q1_dpsi_comp <- run_comparison_test(topmed_hap %>% filter(delta_psi < c1))['bootstrap_p']
r0_q2_dpsi_comp <- run_comparison_test(topmed_hap %>% filter(between(delta_psi, c1, c2)))['bootstrap_p']
r0_q3_dpsi_comp <- run_comparison_test(topmed_hap %>% filter(between(delta_psi, c2, c3)))['bootstrap_p']
r0_q4_dpsi_comp <- run_comparison_test(topmed_hap %>% filter(delta_psi > c3))['bootstrap_p']

#### Delta PSI figure ####
save_plot("fig5S8_delta_psi_bins.svg", width = 6, height = 6)
plot_grid(
  plot_epsilon_plot(r0_q1_dpsi, r0_q1_dpsi_comp) + ggtitle(paste("ΔPSI <", c1)), 
  plot_epsilon_plot(r0_q2_dpsi, r0_q2_dpsi_comp) + ggtitle(paste(c1, "< ΔPSI <", c2)),
  plot_epsilon_plot(r0_q3_dpsi, r0_q3_dpsi_comp) + ggtitle(paste(c2, "< ΔPSI <", c3)),
  plot_epsilon_plot(r0_q4_dpsi, r0_q4_dpsi_comp) + ggtitle(paste("ΔPSI >", c3)),
  nrow = 4
)
dev.off()
#####

# Plot the haplotypes in each group
get_quantile <- function(x){
  ifelse(x < c1, "Q1", 
         ifelse(between(x, c1, c2), "Q2", 
                ifelse(between(x, c2, c3), "Q3", "Q4")))
}

save_plot("fig5S9_vars_per_dpsi_quartile.svg", height = 4, width = 4)
topmed_hap %>%
  select(gene, delta_psi, csnp, CADD) %>%
  mutate(dpsi_quantile = get_quantile(delta_psi)) %>%
  mutate(deleterious = ifelse(CADD > 15, "delterious", "non-deleterious")) %>%
  group_by(deleterious, dpsi_quantile) %>%
  summarize(n = n()) %>%
  ggplot(aes(dpsi_quantile, n, fill = deleterious)) + 
  geom_bar(stat = "identity",position = "dodge") + 
  scale_fill_brewer("CADD Designation", palette = "Set1") +
  ggtitle("Number of haplotypes per delta PSI quantile") +
  xlab("ΔPSI Quartile") +
  gtex_v8_figure_theme()
dev.off()

# Another side note
# Is delta PSI correlated with constraint?
# I'd expect genes with higher delta PSI to be more tolerant to loss-of-function
# variants, and thus have higher LOEUF scores
topmed_hap_genes <- 
  topmed_hap %>%
  select(gene, delta_psi, oe_lof) %>%
  distinct

### Supplemental figure 1 ####
library(ggpubr)
save_plot("fig5S10_dPSI_by_LOEUF_scatter.png", save_fn = png, width = 4, height = 4, units = "in", res = 400)
ggplot(topmed_hap_genes, aes(delta_psi, oe_lof)) + 
  geom_point() + 
  geom_smooth() + 
  gtex_v8_figure_theme() +
  ylab("LOEUF Score") + 
  xlab("ΔPSI") +
  stat_cor(method = "spearman", label.x = .4, cor.coef.name = "rho")
dev.off()


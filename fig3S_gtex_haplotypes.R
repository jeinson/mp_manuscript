# This script is used to perform haplotype enrichment analyses on the GTEx
# data, using more stringent quality filtering and better plotting techniques. 
# 

rm(list = ls())
source("tompen_utility_functions_manuscript.R")
library("stampen") # Can be installed from github
library("ggplot2")
library("magrittr")

# Read in data, using a utility function
gtex_haps <- read_haps("data/fig3_data/gtex_haplotypes/", "secondary_exon_hap_combos_phaser_3_24_22")

# In its current state, how many times do each rare variant show up?
barplot(table(table(gtex_haps$csnp)))

# 1) Filtering ----

# Use SNPs that have read backed phasing or have more than 2 occurences 
csnp_table <- table(gtex_haps$csnp)
non_singletons <- names(csnp_table[csnp_table > 1])
rare_variants <- 
  gtex_haps %>%
  filter(csnp_af_gnomad < .005) %$%
  csnp %>% unique
read_backed_singletons <- 
  gtex_haps %>% 
  filter(csnp_phase_block_snps != ".") %>%
  filter(csnp_cnt == 1) %$% csnp %>% unique

usable_csnps <- intersect(union(non_singletons, read_backed_singletons), rare_variants)

# How many singletons are read backed?
gtex_haps %>%
  filter(csnp_cnt == 1) %$%
  csnp_phase_block_snps %>%
  equals(".") %>%
  table %>% prop.table

# Add a column for the exon ID. 
exon_id_map <- readRDS("data/gtex_v8_exon_id_map.rds")
gtex_haps$top_exon_coord <- exon_id_map[gtex_haps$top_exon_coord]

# Add this so it doesn't break the downstream pipelines
gtex_haps$exp_beta <- gtex_haps$exp_beta_empir

# Addition 3/3/22: Limit to individuals with European ancestry
euro_indvs <- read_tsv("data/gtex_v8_indvs_races.txt") %>%
  deframe
euro_indvs <- names(euro_indvs[euro_indvs == 3])
gtex_haps <- gtex_haps %>% filter(indv %in% euro_indvs)

# F1)
gtex_haps <- gtex_haps %>% filter(csnp %in% usable_csnps)

# F2) Get rid of haplotypes from long genes
haps_per_gene <- table(gtex_haps$gene)

# save_plot("fig3S2_n_haps_per_gene_gtex.svg", width = 5, height = 3)
# boxplot(as.numeric(haps_per_gene), horizontal = T, log = 'x',
#         xlab = "Number of haplotypes per gene in GTEx individuals")
# dev.off()

usable_genes <- flag_outliers(haps_per_gene)
usable_genes <- names(usable_genes)[!usable_genes]
gtex_haps <- gtex_haps %>% filter(gene %in% usable_genes)

# save_plot("fig3S3_n_haps_per_indv_GTEx_V8.svg", width = 4, height = 4)
# hist(table(gtex_haps$indv), breaks = 20, col = "cornflowerblue",
#      main = "", xlab = "N cSNP-sQTL haplotypes per individual")
# box()
# dev.off()

#### Make the 5 summary plots for the supplemenmt ####
## High inclusion allele frequency 
hi_af <- 
  gtex_haps %>%
  select(gene, recom_esnp_af) %>%
  distinct %>%
  ggplot(aes(recom_esnp_af)) +
  geom_histogram(fill = "lightgrey", color = "black") +
  xlab("Frequency of higher exon inclusion allele") +
  ylab("ψQTLs") + 
  gtex_v8_figure_theme() 

## Haplotypes per individual
x <- gtex_haps %>%
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
x <- gtex_haps %>%
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
x <- gtex_haps %>%
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
save_plot("figS3_gtex_haplotype_summaries.svg", width = 10, height = 2)
plot_grid(hi_af, hap_per_indv, hap_p_gene, cv_af, nrow = 1)
dev.off()


# Plot the numbers of haplotypes in a pretty figure
source("scripts/image_xlabel_script.R")
# Note EBImage must be installed from bioconductor for this to work
hap_pic_files <- list.files("data/fig3_data/hap_images/", 
                            full.names = T, pattern = "[01-12].*png")

npoints <- length(hap_pic_files)

### Save the images into a list
hap_pics <- vector(mode = "list", length = npoints)
for(i in 1:npoints){
  hap_pics[[i]] <- EBImage::readImage(hap_pic_files[i])
}

gtex_haps_plt <- gtex_haps
gtex_haps_plt$haplotype <- factor(gtex_haps$haplotype, levels = names(beta_config_sqtl))

collapse_key <- c("abAB", "abAB", "abaB", "abaB", "AbaB", "AbaB", "AbAB", "AbAB")
names(collapse_key) <- names(beta_config_sqtl)

gtex_haps_plt$collapsed_haplotpye <- collapse_key[gtex_haps_plt$haplotype]
gtex_haps_plt$collapsed_haplotpye <- factor(gtex_haps_plt$collapsed_haplotpye, 
                                        levels = c("AbaB", "AbAB", "abAB", "abaB"))

hap_barplot <- 
  ggplot(gtex_haps_plt, aes(collapsed_haplotpye)) + 
  geom_bar(fill = "cornflowerblue") + 
  gtex_v8_figure_theme() +
  theme(plot.margin = unit(c(0.3,0.5,2,0.5), "lines"),
        axis.text.x = element_blank()) 

save_plot("fig3B_gtex_hap_barplot.svg", width = 2.5, height = 1.5)
add_images_as_xlabels(hap_barplot, hap_pics[c(1,3,5,7)]) 
dev.off()

# Make a summary plot of what the number of haplotypes per individual look like
# save_plot("fig3S3_n_haps_per_indv_GTEx_V8.svg", width = 4, height = 4)
# hist(table(gtex_haps$indv), 
#      breaks = 20,
#      main = "",
#      xlab = "N cSNP-sQTL Haplotypes\nper individual", col = "cornflowerblue")
# abline(v = median(table(gtex_haps$indv)), lty = 2, lwd = 1.2)
# box()
# dev.off()

# Run Enrichment on all variants, without splitting by deleterious status. 
bootstrap_test(gtex_haps, exp_beta = "exp_beta_empir")
poibin::ppoibin(sum(gtex_haps$beta), gtex_haps$exp_beta_empir, method = "RNA") * 2

# Check the PSI Z-scores of variants in this set
rv_psi_Z <- readRDS("gtex_lof_in_exons/data/rare_variants_exon_psi_zscore_bytissue.rds") %>%
  bind_rows(.id = "tissue")

sQTLs <- read_tsv("sQTL_v8_anno/cross_tissue_top_sQTLs/top_sQTLs_MAF05.tsv")

gtex_haps_w_psiZ <- 
gtex_haps %>%
  left_join(select(sQTLs, group, tiss), 
            by = c("gene" = "group")) %>%
  left_join(select(rv_psi_Z, indv, tissue, var, exon_psi_Zscore), 
            by = c("csnp" = "var", "tiss" = "tissue", "indv" = "indv")) %>%
  filter(!is.na(exon_psi_Zscore)) %>%
  mutate(deleterious_status = ifelse(CADD > 15, "deleterious", "non-deleterious"))

with(gtex_haps_w_psiZ, tapply(exon_psi_Zscore, deleterious_status, mean, na.rm = T))
t.test(exon_psi_Zscore ~ deleterious_status, gtex_haps_w_psiZ)

# Set haplotypes as factors
gtex_haps$haplotype <- factor(gtex_haps$haplotype, levels = names(beta_config_sqtl))

# Plot 1 - All haplotypes by CADD ----
library(ggplot2)
r0 <- run_enrichment(gtex_haps)
r0_comp <- run_comparison_test(gtex_haps)

# Set to match the TOPMed figure
pval_pos = .023
n_pos = 0.028

save_plot("fig5A_gtex_haplotype_enrichment.svg", width = 5, height = 1.25)
plot_epsilon_plot(r0, pval_pos = pval_pos, n_pos = n_pos,comp_pval = r0_comp['bootstrap_p'])
dev.off()

# What if we look at singletons only?
r1 <- run_enrichment(gtex_haps %>% filter(csnp_cnt == 1))
r1_comp <- run_comparison_test(gtex_haps %>% filter(csnp_cnt == 1))
plot_epsilon_plot(x = r1, comp_pval = r1_comp['bootstrap_p'], pval_pos = .04, n_pos = .045) + 
  ggtitle("Singletons only in GTEx") + ylim(-.05, .05)

# What if we limit to variants on normally highly included exons?
median_psi <- read_tsv("sQTL_v8_anno/multiple_exon_analysis/cross_tissue_top_sQTLs_mean_psi.tsv")
gtex_haps_w_median_psi <- 
  gtex_haps %>%
  left_join(median_psi, by = c('exon_coord' = 'pid'))
r2 <- run_enrichment(gtex_haps_w_median_psi %>% filter(median_psi > .7))
r2_comp <- run_comparison_test(gtex_haps_w_median_psi %>% filter(median_psi > .7))

plot_epsilon_plot(r2, pval_pos = .02, n_pos = .025,comp_pval = r1_comp['bootstrap_p']) + 
  ggtitle("Haplotypes with mean PSI > .7 in GTEx")

# Script for making simple figures that describe the data, which will be used
# for figure 3, second panel. 
# 
# 1/25/22

source("~/myPackages.R")
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance")
source("mp_manuscript/tompen_utility_functions_manuscript.R")

# Figure 1B: High inclusion sQTL allele frequcy distribution

sqtls <- read_tsv("sQTL_v8_anno/cross_tissue_top_sQTLs/top_sQTLs_MAF05.tsv")

library(ggplot2)
sqtl.af.hist <-
  ggplot(sqtls, aes(sQTL_af)) + 
  geom_histogram(fill = "cornflowerblue", 
                 breaks = seq(.05, .95, by = 1/20)) + 
  xlab("Frequency of higher inclusion sQTL allele") +
  ylab("Count") +
  gtex_v8_figure_theme()

save_plot("fig1A_sqtl_af_dist.svg", width = 2.5, height = 1.5)
sqtl.af.hist
dev.off()

# K.S. test for uniformity
ks.test(sqtls$sQTL_af, "punif")
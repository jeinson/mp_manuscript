# Script for making simple figures that describe the data, which will be used
# for figure 3, second panel. 
# 
# 1/25/22

library(tidyverse)
# Set directory to wherever this is saved
source("tompen_utility_functions_manuscript.R")

# Figure 1B: High inclusion sQTL allele frequcy distribution

# See https://www.biorxiv.org/content/10.1101/2022.12.05.518915v1 for details
# on this file
sqtls <- read_tsv("data/fig3_data/top_sQTLs_MAF05.tsv")

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
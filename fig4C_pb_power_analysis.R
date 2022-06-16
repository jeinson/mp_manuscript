# Plot the results from the power analysis
# Jonah Einson
# 
# 2/26/22

rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/analysis/bootstrap_poisson_binomial")
source("~/myPackages.R")
source("../../mp_manuscript/tompen_utility_functions_manuscript.R")
library(tidyverse)

res_fp <- list.files("power_results", pattern = "fine", full.names = T)

res <- map_dfr(res_fp, read_tsv)
res$`Beta Parameters` <- factor(with(res, paste0("⍺ = ", alpha, ", β = ", beta)))

# To make sure we only use the results from the right betas
params <- data.frame(
  alpha = c(1, 2, .3, 50, 1.387),
  beta = c(1, 6, .3, 50, 0.954)
)
res <- inner_join(res, params)

fig2c <- 
  ggplot(res, aes(sd_from_mean, BS, color = `Beta Parameters`)) + 
  #geom_point() + 
  geom_line() + 
  xlab("Standard deviations from the mean") +
  ylab("Power") +
  gtex_v8_figure_theme() +
  theme(panel.grid.major = element_line(color = "grey", linetype = 2, size = .5)) + 
  geom_hline(yintercept = 0, lty = 1, col = "grey") + 
  geom_vline(xintercept = 0, lty = 1, col = "darkgrey")

# Save plot to disk as svg
base_path = "/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/mp_manuscript"
svg(paste0(base_path, "/figures/fig2C.svg"), width = 4, height = 2)
fig2c
dev.off()

fig2d <- 
  ggplot(res, aes(sd_from_mean, mean_epsilon, color = `Beta Parameters`)) + 
  #geom_point() + 
  geom_line() + 
  xlab("Standard deviations from the mean") +
  ylab("ε") +
  gtex_v8_figure_theme() +
  theme(panel.grid.major = element_line(color = "grey", linetype = 2, size = .5)) + 
  geom_hline(yintercept = 0, lty = 1, col = "grey") + 
  geom_vline(xintercept = 0, lty = 1, col = "darkgrey")

svg(paste0(base_path, "/figures/fig2D.svg"), width = 4, height = 2, )
fig2d
dev.off()


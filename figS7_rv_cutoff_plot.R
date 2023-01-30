# Figure S7: Plot enrichment figures with subsequently increasing allele
# frequency cutoffs for TOPMed data. 

rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/")
source("~/myPackages.R")
library(here)
source("mp_manuscript/tompen_utility_functions_manuscript.R")
library(tompen)
library(ggplot2)
library(cowplot)

# Read in results from supplement_fig_minor_allele_cutoff_on_hap_depletion.R
r2 <- readRDS("topmed/cutoff_results/r2.rds")
r3 <- readRDS("topmed/cutoff_results/r3.rds")
r4 <- readRDS("topmed/cutoff_results/r4.rds")
r5 <- readRDS("topmed/cutoff_results/r5.rds")

r2_comp <- readRDS("topmed/cutoff_results/r2_comp.rds")['bootstrap_p']
r3_comp <- readRDS("topmed/cutoff_results/r3_comp.rds")['bootstrap_p']
r4_comp <- readRDS("topmed/cutoff_results/r4_comp.rds")['bootstrap_p']
r5_comp <- readRDS("topmed/cutoff_results/r5_comp.rds")['bootstrap_p']

n_pos = .02
pval_pos = .015
plot_grid(
  plot_epsilon_plot(r2, comp_pval = r2_comp, pval_pos = pval_pos, n_pos = n_pos),
  plot_epsilon_plot(r3, comp_pval = r3_comp, pval_pos = pval_pos, n_pos = n_pos), 
  plot_epsilon_plot(r4, comp_pval = r4_comp, pval_pos = pval_pos, n_pos = n_pos), 
  plot_epsilon_plot(r5, comp_pval = r5_comp, pval_pos = pval_pos, n_pos = n_pos), 
  ncol = 1
)

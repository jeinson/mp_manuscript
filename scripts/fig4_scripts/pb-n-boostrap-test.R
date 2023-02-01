# This script is used to benchmark the number of iterations needed to accurately
# estimate the poisson-binomial distribution using the bootstrap method. 
# 
# Jonah Einson
# 10/12/21

library(poibin)
source("scripts/fig4_scripts/pb-functions.R")

set.seed(1234567890)

# Run the code
n = 1000
test_alpha = 1.378
test_beta = 0.954
dat <- create_sample_dataset(n, test_alpha, test_beta)

Fx_DFT <- ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "DFT-CF")

n_boots <- c(100, 1000, 2500, 5000, 10000, 15000)
res_list <- list()

# Test this with different numbers of bootstraps
for(boot in n_boots){
  message(paste("Checking with", boot, "bootstrapped resamplings\n"))
  res_list[[as.character(boot)]] <- 
    replicate(50, ppoibin_bootstrap(kk = dat$kk, pp = dat$pp, B = boot))
}

# Reformat the results for plotting
library(tidyverse)
results_df <- 
  res_list %>%
  map(as.data.frame) %>%
  bind_rows(.id = "N") %>%
  rename("X" = ".x[[i]]")

results_df$diff <- results_df$X - Fx_DFT
results_df$N <- factor(results_df$N, levels = n_boots)

plt <- ggplot(results_df, aes(N, diff)) + 
  geom_violin() + 
  geom_hline(yintercept = 0) + 
  xlab("N bootstraps") + 
  ylab("Error") + 
  theme_bw()

png("figures/n_bootstraps_analysis.png", width = 4, height = 4, res = 300, units = "in")
plt
dev.off()

# Test this again in relation on a simple binomial distributed RV
n = 100
pp_bin = rep(.5, n)
kk_bin = c(rep(1, n/2), rep(0, n/2))
x <- replicate(100, 1 - ppoibin_bootstrap(kk_bin, pp_bin, B = 1000))

boxplot(x, horizontal = T)
abline(v = pbinom(q = sum(kk_bin), length(pp_bin), .5), col = "red")


# This script reads in the results from pb-test-benchmark.R and generates
# tables and compiles some figures. 
# 
# Jonah Einson
# jeinson@nygenome.org
# 10/14/21

setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/analysis/bootstrap_poisson_binomial")
source("../../mp_manuscript/tompen_utility_functions_manuscript.R")
source("~/myPackages.R")
library(ggplot2)

# Read in results
fp <- list.files("results", full.names = T, pattern = ".tsv")

params <- data.frame(
  alpha = c(1, 2, .3, 50, 1.387),
  beta = c(1, 6, .3, 50, 0.954)
)
params$beta_dist_label <- factor(with(params, paste0("ɑ = ", alpha, ", β = ", beta)))

results <- lapply(fp, function(x){
  y = unlist(str_split(x, "_"))[c(5,7,9)]
  df = read_tsv(x)
  df = colMeans(abs(df - df$DFT))
  df['alpha'] = as.numeric(y[2])
  df$beta = as.numeric(str_remove(y[3], ".tsv"))
  df$N = as.numeric(y[1])
  
  df
})

results <- results %>% map(as.data.frame) %>% bind_rows()

# Make some plots
results$beta_dist_label <- factor(with(results, paste0("ɑ = ", alpha, ", β = ", beta)))
results <- results %>% arrange(beta_dist_label, N)
results$logN <- round(log10(results$N), digits = 1)
results$logN <- fct_inorder(as.character(results$logN))

# Make sure we're only using the ones we need. 
results <- results %>% filter(beta_dist_label %in% params$beta_dist_label)

df <- gather(results, key = "Estimation Method", value = "Mean Error", 
             c("RNA", "Na", "PA", "BS", "Naive"))
df$`Estimation Method` <- factor(df$`Estimation Method`, 
                                 levels = c("BS", "PA", "Na", "RNA", "Naive"))

plt <- 
  ggplot(df, aes(logN, `Mean Error`, col = `Estimation Method`, group = `Estimation Method`)) + 
  geom_point() + 
  stat_summary(fun = sum, geom = "line", lwd = 1) + 
  facet_wrap(facet = "beta_dist_label") + 
  theme_bw() + 
  theme(legend.position = c(0.9, 0.15)) + 
  xlab(bquote(log[10]~"("~Simulated~samples~")")) +
  scale_color_brewer(
    name = "Estimation Method", 
    labels = c("Bootstrap", 
               "Poisson Approximation", 
               "Normal Approximation", 
               "Regularized Normal Approximation", 
               "Naive approach"), palette = "Dark2") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=.5))

#png("figures/mean_error_by_sample_size.png", width = 9, height = 4, res = 300, units = "in")
save_plot("fig4S1_poisson_binomial_benchmark.svg", width = 6, height = 4)
plt
dev.off()

# Also make a plot showing how the Poisson-binomial correlates with FFT based
# method for calculating the PB. 
res_raw <- read_tsv("results/results_N_1000_n_1000_alpha_1.387_beta_0.954.tsv")
res_pval <- data.frame(
  BS_p = -log10(sapply(res_raw$BS, function(x) min(x, 1-x)) * 2), 
  DFT_p = -log10(sapply(res_raw$DFT, function(x) min(x, 1-x)) * 2)
  )

# Get rid of the infitinte rows
res_pval <- res_pval[is.finite(rowSums(res_pval)),]

save_plot("fig4C_poisson_binomial_DFT_raw_comparison.svg", width = 2, height = 2)
ggplot(res_pval, aes(BS_p, DFT_p)) + 
  geom_point(shape = 1) + 
  xlab("-log10(Bootstrap test p-value)") + 
  ylab("-log10(Poisson-binomial p-value") + 
  geom_smooth(method = "lm", col = "red", lwd = .5) +
  gtex_v8_figure_theme()
dev.off()  

# Get the correlation and stuff for the figure

lin_mod <- lm(DFT_p ~ BS_p, data = res_pval)
summary(lin_mod)

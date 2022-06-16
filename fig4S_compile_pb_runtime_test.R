# Compile the results from the runtime benchmark. 
# 
# Jonah Einson
# 2/14/22 (HVD)

setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/analysis/bootstrap_poisson_binomial")
source("~/myPackages.R")

results_list <- list()
# N = c(100, 200, 500, 1000, 2000, 5000, 10000, 15000, 25000)#, 50000)
N = floor(10 ^ seq(1, 4.5, by = .5))
for(n in N){
  results_list[[as.character(n)]] <- 
    readRDS(paste0("runtime_results/", n, "_replicate_results.rds"))
}

runtime_df <- 
  results_list %>%
  map(as.data.frame) %>%
  bind_rows(.id = "N") %>%
  mutate(N = as.numeric(N)) %>%
  mutate(logN = round(log10(N), 1)) %>%
  mutate(logN = factor(logN)) %>%
  filter(expr != "RF")

library(scales)

base_path = "/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/mp_manuscript/"
# svg(paste0(base_path, "figures/fig2S_pb_runtime.svg"), width = 8, height = 2)
ggplot(runtime_df, aes(logN, time * 10e-9, color = expr)) +
  geom_boxplot() + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("Time (seconds)") +
  xlab(expression(log[10]~N~simulated~datapoints)) +
  scale_color_brewer(
    name = "Estimation Method", 
    labels = c("Bootstrap Method",
               "Poisson Approximation",
               "Normal Approximation",
               "Regularized Normal Approximation",
               "Naive approach"),
    palette = "Set2") +
  theme_linedraw()
# dev.off()


######
# This script calculates the PDF of various P-B distributed random samples, using
# different methods, in order to evaluate the accuracy under various conditions
#
# Jonah Einson
# 10/12/2021

setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/analysis/bootstrap_poisson_binomial")
library(poibin)
source("~/myPackages.R")
source("pb-functions.R")

# Parse arguments from the command line
if(interactive()){
  replicates = 1000
  test_alpha =2
  test_beta = 6
  n = 100
} else {
  args <- commandArgs(TRUE)
  replicates = as.numeric(args[1])
  test_alpha = as.numeric(args[2])
  test_beta = as.numeric(args[3])
  n = as.numeric(args[4])
}

# Set a seed for re-doing this later
set.seed(replicates + test_alpha + test_beta)

# Calculate PPs using different PB approximation methods
N = replicates
# ns = c(10, 20, 50, 100, 200, 500, 
#        1000, 2000, 5000, 10000, 150000)

result_list <- list()

# for(n in ns){
message(paste("Samples drawn:", n))
out_mat <- matrix(nrow = N, ncol = 7)

for(i in 1:N){
  if(i %% 10 == 0) message(paste0(i, "..."), appendLF = F)
  dat <- create_sample_dataset(n, test_alpha, test_beta)
  Fx_DFT <- ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "DFT-CF")
  Fx_RF <- ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "RF")
  Fx_RNA <- ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "RNA")
  Fx_NA <- ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "NA")
  Fx_PA <- ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "PA")
  
  # My package that simulates the pois-binom distribution
  Fx_BS <- ppoibin_bootstrap(kk = dat$kk, pp = dat$pp, B = 5000)
  
  # Do this again with the naive approach
  Fx_naieve <- ppoibin_naieve(kk = dat$kk, pp = dat$pp)
  
  results <- c("DFT" = Fx_DFT, "RF" = Fx_RF, "RNA" = Fx_RNA, 
               "Na" = Fx_NA, "PA" = Fx_PA,
               "BS" = Fx_BS, "Naive" = Fx_naieve)
  
  out_mat[i, ] <- results
}

out_mat <- as.data.frame(out_mat)
colnames(out_mat) <- names(results)
result_list[[as.character(n)]] <- colMeans(abs(out_mat - out_mat$DFT))
message("\n")
# }

out_df <- dplyr::bind_rows(result_list, .id = "N")

# Add other info back to the dataframe. I'll thank myself later
out_df$alpha <- test_alpha
out_df$beta <- test_beta

out_name <- paste0("results/results_N_", replicates, 
                   "_n_", n, 
                   "_alpha_", test_alpha, 
                   "_beta_", test_beta, ".tsv")

# readr::write_tsv(out_df, out_name)

# Actually save the whole thing. easier to use downstream
readr::write_tsv(out_mat, out_name)

# Make the sample sheet. This will be commented out after it's run once. 
# params <- data.frame(
#   alpha = c(1, 2, .3, 50, 1.387),
#   beta = c(1, 6, .3, 50, 0.954)
# )
# 
# replicates = 1000
# # sample_size = c(10,20,50,100,200,500,1000,2000,5000,10000,15000)
# sample_size = floor(10 ^ seq(1, 4.5, by = .5))
# 
# merge(merge(replicates, params), sample_size) %>%
#   arrange(alpha, beta) %>%
#   set_colnames(c("replicates", "alpha", "beta", "n")) ->
#   sample_sheet
# 
# write_tsv(sample_sheet, "sample_sheets/poisson_binomial_bootstrap_sample_sheet.txt")

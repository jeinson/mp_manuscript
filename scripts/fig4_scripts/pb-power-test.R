# This script runs some power analyses of the test on different versions

library(poibin)
library(purrr)
library(readr)
source("scripts/fig4_scripts/pb-functions.R")

# Parse arguments from the command line
if(interactive()){
  test_alpha = .3
  test_beta = .3
  sd_from_mean = -.4
  replicates = 100
} else {
  args <- commandArgs(TRUE)
  test_alpha = as.numeric(args[1])
  test_beta = as.numeric(args[2])
  sd_from_mean = as.numeric(args[3])
  replicates = 1000
}

n = 1000

# Set a seed for re-doing this later
set.seed(replicates + test_alpha + test_beta)

# Get the mean and variance of the distribution used to simulate data
beta_mean <- test_alpha / (test_alpha + test_beta)
beta_sd <- sqrt(
  (test_alpha * test_beta) / 
    ((test_alpha + test_beta)^2 * (test_alpha + test_beta + 1))
  )

# Calculate PPs using different PB approximation methods
N = replicates
out_mat <- matrix(nrow = N, ncol = 3)
for(i in 1:N){
  if(i %% 10 == 0) message(paste0(i, "..."), appendLF = F)
  dat <- create_sample_dataset(n, test_alpha, test_beta)
  
  # Modify the data so the rate of success of kk is simulated_p_hp
  simulated_p_hp <- beta_mean + (sd_from_mean * beta_sd)
  dat$kk <- as.numeric(rbernoulli(length(dat$kk), p = simulated_p_hp))
  
  # Now run the test with the pb-bootstrap method
  # My package that simulates the pois-binom distribution
  Fx_BS <- ppoibin_bootstrap(kk = dat$kk, pp = dat$pp, B = 5000)
  
  # Do this again with the naive approach
  Fx_naieve <- ppoibin_naieve(kk = dat$kk, pp = dat$pp)
  
  # Also save the epsilon score
  epsilon <- with(dat, mean(kk - pp))
  
  results <- c("BS" = Fx_BS, "Naive" = Fx_naieve)
  
  # Convert to a p-value
  results <- 2 * if(results['BS'] < .5){
    results
    } else {
      1-results
    }
  
  # Add the epsilon result in
  results <- c(results, "epsilon" = epsilon)
  out_mat[i, ] <- results
}

out_mat <- as.data.frame(out_mat)
colnames(out_mat) <- names(results)

# Calculate power
power_res <- colSums(out_mat[,1:2] < .05) / nrow(out_mat)

# Mean epsilon score across all the replicates
mean_epsilon <- mean(out_mat$epsilon)

# Save the output
out_df <- cbind(data.frame(
  "alpha" = test_alpha, 
  "beta" = test_beta, 
  "sd_from_mean" = sd_from_mean), 
  t(data.frame(power_res)), 
  "mean_epsilon" = mean_epsilon
)

out_nm <- paste0("power_results/power_res_", test_alpha, "_", test_beta, "_", sd_from_mean, "_fine.txt")
write_tsv(out_df, out_nm)


# Generate a sample sheet (will comment out later)
# params <- data.frame(
#   alpha = c(1, 2, .3, 50, 1.387),
#   beta = c(1, 6, .3, 50, 0.954)
# )
# sd_from_mean <- round(seq(-.1, .1, by = .01), 2)
# sample_sheet <- merge(params, as.data.frame(sd_from_mean)) %>%
#   arrange(alpha, beta)
# write_tsv(sample_sheet, "sample_sheets/pb_power_analysis_sample_sheet_fine.tsv",
#           col_names = F)


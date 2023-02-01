# This script looks at how the P-B test handles simulated hqplotype data from
# the tompen package. This would be very interesting to test for. 
# 
# Jonah Einson
# 11/23/21

library(stampen)
library(poibin)
library(magrittr)

# simulated_haps <- simulate_haplotype_counts(100, 500)
# simulated_haps <- characterize_haplotypes(simulated_haps, beta_config = beta_config_sqtl)

simulated_haps <- simulate_null_haplotypes(100, 500)
simulated_haps <- characterize_haplotypes(simulated_haps, beta_config_sqtl)

# Check AF and exp_beta correlation
simulated_haps %$% 
  plot(sqtl_af, exp_beta, 
       xlab = "Simulated high-inclusion sQTL allele frequency", 
       ylab = "E[β]")


# Now test the various methods, and use the different options for expectations
# of beta. 

# DFT
Fx_DFT_dip <- ppoibin(kk = sum(simulated_haps$beta),
                     pp = simulated_haps$exp_beta, method = "DFT-CF")
Fx_DFT_hap <- ppoibin(kk = sum(simulated_haps$beta),
                      pp = simulated_haps$sqtl_af, method = "DFT-CF")

# RF Method
Fx_RF_dip <- ppoibin(kk = sum(simulated_haps$beta),
                      pp = simulated_haps$exp_beta, method = "RF")
Fx_RF_hap <- ppoibin(kk = sum(simulated_haps$beta),
                      pp = simulated_haps$exp_beta, method = "RF")

# How do the allele frequencies vs. the expectations of beta differ?
plot(simulated_haps$sqtl_af, simulated_haps$exp_beta)

# Repeat this process a few times to check if they are generally the same thing
reps = 200
set.seed(1124215)
# Do this in parallel to speed things up
library(doParallel)
registerDoParallel(8)

out_df <- 
  foreach(icount(reps), .combine=rbind) %dopar% {
    # Generate fake data
    simulated_haps <- simulate_null_haplotypes(100, 500)
    simulated_haps <- characterize_haplotypes(simulated_haps, beta_config = beta_config_sqtl)
    
    # Calculate 
    Fx_DFT_dip <- ppoibin(kk = sum(simulated_haps$beta),
                          pp = simulated_haps$exp_beta, method = "DFT-CF")
    Fx_DFT_hap <- ppoibin(kk = sum(simulated_haps$beta),
                          pp = simulated_haps$sqtl_af, method = "DFT-CF")
    
    # Save results 
    c(dip = Fx_DFT_dip, hap = Fx_DFT_hap, af_df = sd(simulated_haps$sqtl_af))
  }

out_df <- as.data.frame(out_df)
mod = lm(dip~hap, data = out_df)

plot(out_df$dip, out_df$hap, 
     xlab = "p(X = x) using M0", 
     ylab = "p(X = x) using M1")
abline(a = 0, b = 1, col = "red", lty = 3, lwd = 5)
abline(mod, col = "blue")

# Is the SD of the probability vector correlated with the difference?
out_df$diff <- out_df$dip - out_df$hap

out_df %$%
  plot(af_df, diff, 
       xlab = "SD of the simulated allele frequencies", 
       ylab = "Difference p(x) using both approaches")

# dat <- create_sample_dataset(n, test_alpha, test_beta)
# Fx_DFT <- ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "DFT-CF")
# Fx_RF <- ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "RF")
# Fx_RNA <- ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "RNA")
# Fx_NA <- ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "NA")
# Fx_PA <- ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "PA")
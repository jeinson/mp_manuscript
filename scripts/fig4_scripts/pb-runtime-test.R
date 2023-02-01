#!/nfs/sw/R/R-4.0.0/bin/Rscript
# This script is used to compute the runtime of the different methods
# for calculating the PDF of the poisson-binomial distribution. 
# 
# Jonah Einson
# 10/12/21

library(poibin)
library(microbenchmark)

# User Defined Functions
create_sample_dataset <- function(N, a, b){ 
  # N: Sample size
  # a: alpha hyper-perameter for the Beta distribution
  # b: beta hyper-parameter for the Beta distribution
  
  pp <- rbeta(N, a, b)
  kk <- sapply(pp, function(x) sample(c(1,0), 1, prob = c(x, 1-x)))
  
  list("pp" = pp, "kk" = kk)
}

ppoibin_bootstrap <- function(kk, pp, B = 1000) {
  # kk: A vector of binary responses
  # pp: A vector of p(x) for each binary response
  # B: The number of bootstraps (1000 by default)
  epsilon <- mean(kk - pp)
  p_b = rep(0, B)
  for (b in 1:B) {
    ix <- sample(1:length(kk), length(kk), replace = T)
    p_b[b] <- mean(kk[ix] - pp[ix])
  }
  p_b <- sort(p_b)
  sum((p_b > 0) / B)
}

path_fn <- function(n_reps){
  paste0('/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/analysis/bootstrap_poisson_binomial/runtime_results/',
         n_reps, "_replicate_results.rds")
}

# run the code
# Start with 100
#N = c(100, 200, 500, 1000, 2000, 5000, 10000, 15000, 25000, 50000)
N = floor(10 ^ seq(1, 5, by = .5))
for(n in N){
  message(n)
  test_alpha = 1
  test_beta = 1
  dat <- create_sample_dataset(n, test_alpha, test_beta)
  
  mbm <- microbenchmark(
    "DFT-CF" = ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "DFT-CF"),
    "RF"  = ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "RF"),
    "RNA" = ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "RNA"),
    "NA"  = ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "NA"),
    "PA"  = ppoibin(kk = sum(dat$kk), pp = dat$pp, method = "PA"),
    "BS" =  ppoibin_bootstrap(kk = dat$kk, pp = dat$pp, B = 5000), 
    times = 10
  )
  
  readr::write_rds(mbm, path_fn(n))
  
  # plt_title <- paste("Runtime of P-B methods with", n, "samples")
  # plot(mbm, main = plt_title)
  
}

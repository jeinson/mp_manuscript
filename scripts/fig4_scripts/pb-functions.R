# Define two function used to test the Poisson-binomial bootstrap method

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
  sum(p_b > 0) / B
}

# A "Naieve" approach, where we just take the average of the probability 
# vector and use that as the parameter in a binomial distribution
ppoibin_naieve <- function(kk, pp){
  ## kk: A vector of binary responses
  ## pp: A vector of p(x) for each binary response
  
  p <- mean(pp)
  q <- sum(kk)
  n <- length(kk)
  
  # Just run pbinom
  pbinom(q = q, size = n, prob = p)
}

library(matlib)
library(data.table)
library(Matrix)

generate_data_H_0 <- function(n, p, distr_eps) {
  eps <- calculate_eps(n, p, distr_eps)
  B <- matrix(runif(p^2, -1, 1), nrow = p)
  diag(B) <- 1
  X <- t(solve(B) %*% eps)
  return(X)
}

calculate_eps <- function(n, p, distr_eps) {
  if (distr_eps == "gamma") {
    shapes <- runif(p, 1, 3)
    rates <- runif(p, 1, 5)
    eps <- t(mapply(function(n, shape, rate) rgamma(n, shape, rate) - shape / rate, n, shapes, rates))
  } else if (distr_eps == "beta") {
    alphas <- runif(p, 0.5, 2)
    betas <- runif(p, 2, 10)
    eps <- t(mapply(function(n, beta, alpha) rbeta(n, alpha, beta) - alpha/(alpha+beta), n, alphas,  betas))
  } else if (distr_eps == "overlapping_gaussian") {
    sds_1 <- runif(p, 0.5, 1)
    sds_2 <- runif(p, 1, 3)
    mean_diffs <- runif(p, 3, 4)
    eps <- t(mapply(sample_overlapping_gaussian, n, sds_1, sds_2, mean_diffs))
  } else if (distr_eps == "gaussian") {
    sds <- runif(p, 0.5, 2)
    eps <- t(mapply(function(n, sd) rnorm(n, 0, sd), n, sds))
  } else {
    stop("invalid parameter for distr_eps")
  }
  return(eps)
}

sample_overlapping_gaussian <- function(n, sd_1, sd_2, mean_diff) {
  # Sample so that order is random
  return(sample(c(rnorm(n/2, 0, sd_1), rnorm(n/2, mean_diff, sd_2)) - 1/2 * mean_diff))
}
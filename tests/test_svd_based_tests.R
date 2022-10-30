# imports
source("goodness-of-fit-tests\\svd_based_tests.R")
library(testthat)

test_that("M is estimated correctly", {
    # Constant data
    n <- 5
    p <- 2
    X <- matrix(1, n, p)
    M_hat <- estimate_M(X, p)
    expect_equal(M_hat, matrix(1, p+1, p+1))

    # Gaussian data, test s_12
    n <- 1000
    p <- 2
    X <- matrix(rnorm(p*n), n, p)
    M_hat <- estimate_M(X, p)
    expected_s12 <- mean(X[,1]*X[,2])
    tolerance <- 1e-10
    expect_true(abs(M_hat[1,2] - expected_s12) < tolerance)
})

test_that("W is estimated correctly", {
    # Constant data, p = 2
    n <- 500
    p <- 2
    X <- matrix(1, n, p)
    W_hat <- estimate_W(X, p)
    expect_equal(W_hat, matrix(0, (p+1)^2, (p+1)^2))

    # Constant data, p = 3
    n <- 500
    p <- 3
    X <- matrix(1, n, p)
    W_hat <- estimate_W(X, p)
    expect_equal(W_hat, matrix(0, 4 * 6, 4 * 6))
})

test_that("Rank is estimated correctly", {
    # Generate data such that M has lower rank
    n <- 2
    p <- 2
    X <- matrix(1, n, p)
    X[, 2] <- c(1, 0)

    rk <- test_KP(X, p, p)[["TSTAT"]]
    tolerance <- 1e-10
    expect_true(abs(rk - 0) < tolerance)
})
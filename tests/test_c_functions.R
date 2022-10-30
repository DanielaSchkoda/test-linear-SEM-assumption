library(Rcpp)
library(testthat)
source("goodness-of-fit-tests\\minor_based_tests.R")
sourceCpp("goodness-of-fit-tests\\minor_based_tests_c_optimization.cpp")

test_that("polynomial estimates correct", {
    # R function
    poly = list(c(s12=1, s34=1, coef=1), c(s13=1, s24=1, coef=-1), c(t123=1, t234=2, coef=1))
    est_p1 = create_estimator(poly)
    L = list(c(2,4,6,8), c(3,5,7,1), c(3,4,2,5))
    # expected: 67196
    expect_equal(est_p1(L), 67196)

    # C++ function
    poly_c = polynomial_to_c_format(poly)
    L_c = do.call("rbind", L)
    expect_equal(estimate_polynomial_c(poly_c, L_c), 67196)
})

test_that("H is calculated correctly", {
    poly1 = list(c(s12=1, s34=1, coef=1), c(s13=1, s24=1, coef=-1), c(t123=1, t234=2, coef=1))
    poly2 = list(c(s13=1, s24=1, coef=1), c(s14=1, s23=1, coef=-3), c(t124=1, t234=1, s11=1, coef=2))
    polynomials = list(poly1, poly2)
    X = matrix(sample(-5:5, 4*6, replace = TRUE), 6, 4)
    indices = t(sapply(1:8, function(i) sample(1:6, 4, replace = FALSE)))

    H_R = H(X, indices, lapply(polynomials, create_estimator))
    H_C = H_c(X, indices, lapply(polynomials, polynomial_to_c_format))
    tolerance <- 1e-10
    expect_true(all(abs(H_R - H_C) < tolerance))
})
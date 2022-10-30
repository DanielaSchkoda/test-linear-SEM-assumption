library(data.table)
library(matlib)
library(gtools)
library(Matrix)
library(expm)
library(wrapr)
library(boot)

generate_pairs <- function(p) {
    pairs <- combinations(p, 2, 1:p, repeats.allowed = T)
    pairs <- lapply(split(pairs, seq(nrow(pairs))), function(ind) unlist(ind, use.names = FALSE))
    return(pairs)
}

estimate_moment <- function(X, ind) {
    X_components_multiplied <- apply(X[, ind], 1, prod)
    return(mean(X_components_multiplied))
}

estimate_M <- function(X, p) {
    column_indices <- generate_pairs(p)

    M <- matrix(, nrow = 1 + p, ncol = length(column_indices))
    for (i in 1:nrow(M)) {
        for (j in 1:length((column_indices))) {
            if (i == 1) {
                M[i, j] <- estimate_moment(X, column_indices[[j]])
            } else {
                M[i, j] <- estimate_moment(X, append(column_indices[[j]], i - 1))
            }
        }
    }
    return(M)
}

estimate_W <- function(X, p) {
    column_indices_M <- generate_pairs(p)
    indices_M_by_col <- lapply(column_indices_M, function(col_ind) c(list(col_ind), lapply(1:p, function(row_ind) append(col_ind, row_ind))))
    indices_W <- unlist(indices_M_by_col, recursive = FALSE)
    W <- matrix(, length(indices_W), length(indices_W))
    for (i in 1:nrow(W)) {
        for (j in 1:ncol(W)) {
            W[i, j] <- estimate_moment(X, append(indices_W[[i]], indices_W[[j]])) - estimate_moment(X, indices_W[[i]]) * estimate_moment(X, indices_W[[j]])
        }
    }
    return(W)
}

calculate_A_perp_B_perp <- function(M, r) {
    # svd returns S, U, V s.t. M = U * S * V"
    svd_decomp <- svd(M, nu = nrow(M), nv = ncol(M))
    U_22 <- as.matrix(svd_decomp[[2]][(r + 1):nrow(M), (r + 1):nrow(M)])
    V_22 <- as.matrix(svd_decomp[[3]][(r + 1):ncol(M), (r + 1):ncol(M)])

    A_perp <- svd_decomp[[2]][, (r + 1):nrow(M)] %*% solve(U_22) %*% sqrtm(U_22 %*% t(U_22))
    B_perp <- sqrtm(V_22 %*% t(V_22)) %*% solve(t(V_22)) %*% t(svd_decomp[[3]][, (r + 1):ncol(M)])
    return(list("A_perp"=A_perp, "B_perp"=B_perp))
}

##################
## KP rank test ##
##################

test_KP <- function(X, p, r, scaling=FALSE) {
    n <- nrow(X)
    M <- estimate_M(X, p)

    # Scaling
    if (scaling) {
      G <- matrix(0, nrow = nrow(M), ncol = nrow(M))
      diag(G) <- c(1, 1/sqrt(apply(X^2, 2, mean)))
      column_indices <- generate_pairs(p)
      H <- matrix(0, nrow = ncol(M), ncol = ncol(M))
      diag(H) <- sapply(column_indices, function(ind) 1/sqrt(mean(X[,ind[1]]^2)) * 1/sqrt(mean(X[,ind[2]]^2)))
    } else {
      G <- diag(nrow = nrow(M))
      H <- diag(nrow = ncol(M))
    }

    Theta <- G %*% M %*% t(H)

    calculate_A_perp_B_perp(Theta, r) %.>% to(
        A_perp <- A_perp,
        B_perp <- B_perp
    )
    Lambda_q <- matrix(t(A_perp) %*% Theta %*% t(B_perp), ncol = 1, byrow=FALSE)

    W <- (H %x% G) %*% estimate_W(X, p) %*% t(H %x% G)
    Omega <- (B_perp %x% t(A_perp)) %*% W %*% t(B_perp %x% t(A_perp))

    # That's wrong in paper it needs to be n not 1/n. See previous step in paper or https://github.com/miabrahams/UsefulMatlabCode/blob/master/KP_RankTest.m
    rk <- n * (t(Lambda_q) %*% solve(Omega) %*% Lambda_q)[1, 1]
    df <- (nrow(M) - r) * (ncol(M) - r)
    pval <- 1 - pchisq(q = rk, df = df)

    return(list("PVAL" = pval, "TSTAT" = rk, "estimated_M" = M))
}

#######################
# Test with bootstrap #
#######################

calculate_test_stat_svd <- function(data, ind, rank, M_hat, U_2, V_2) {
    n <- nrow(data)
    p <- ncol(data)
    bootstrap_sample <- data[ind, ]
    M_bootstrap <- estimate_M(bootstrap_sample, p)
    mathcal_M <- sqrt(n) * (M_hat - M_bootstrap)
    phi_r <- sum(svd(t(U_2) %*% mathcal_M %*% V_2)[[1]][1:(nrow(M_hat) - rank)]^2)
    return(phi_r)
}

test_svd_bootstrap <- function(X, r, E=1000) {
    n <- nrow(X)
    p <- ncol(X)
    M_hat <- estimate_M(X, p)

    # function svd returns S, U, V s.t. M = U * S * V"
    svd_decomp <- svd(M_hat, nu = nrow(M_hat), nv = ncol(M_hat))
    phi_r <- sum(svd_decomp[[1]][(r + 1):nrow(M_hat)]^2)

    U_2 <- as.matrix(svd_decomp[[2]][, (r + 1):nrow(M_hat)])
    V_2 <- as.matrix(svd_decomp[[3]][, (r + 1):ncol(M_hat)])

    # bootstrapping
    bootstrap_samples <- boot(data = X, statistic = calculate_test_stat_svd, R = E, rank = r, M_hat = M_hat, U_2 = U_2, V_2 = V_2)$t

    pval <- (1+length(which(c(bootstrap_samples) >= n * phi_r)))/(1+E)

    return(list("PVAL" = pval, "TSTAT" = phi_r))
}

# SVD bootstrap test that takes possibility of lower rank into account
test_svd_bootstrap_lower_rank <- function(X, p, r, E=1000) {
    n <- nrow(X)
    M_hat <- estimate_M(X, p)
    alpha <- 0.05
    beta <- alpha/10

    # find r_0
    r_0 <- r+1 
    for (r_prime in 1:r) {
        # If one of the tests accepts, set r_0 = r
        res <- test_KP(X, p, r_prime,scaling=TRUE)
        if(res[["PVAL"]] >= beta) {
            r_0 <- r_prime
            break
        }
    }

    if (r_0 == r+1) {
        # reject already and return PVAL, TSTAT of last KP iteration
        return(list("PVAL" = res[["PVAL"]], "TSTAT" = res[["TSTAT"]], "r_0" = r_0))
    }

    # svd returns S, U, V s.t. M = U * S * V'
    svd_decomp <- svd(M_hat, nu = nrow(M_hat), nv = ncol(M_hat))
    phi_r <- sum(svd_decomp[[1]][(r + 1):nrow(M_hat)]^2)

    U_2 <- as.matrix(svd_decomp[[2]][, (r + 1):nrow(M_hat)])
    V_2 <- as.matrix(svd_decomp[[3]][, (r + 1):ncol(M_hat)])

    # bootstrapping
    bootstrap_samples <- boot(data = X, statistic = calculate_test_stat_svd, R = E, rank = r, M_hat = M_hat, U_2 = U_2, V_2 = V_2)$t

    pval <- (1+length(which(c(bootstrap_samples) >= n * phi_r)))/(1+E)

    return(list("PVAL" = pval, "TSTAT" = n * phi_r, "r_0" = r_0))
}
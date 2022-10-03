# The code is adapted from Nils Sturma, https://github.com/NilsSturma/TestGGM

# imports
library(gtools)
library(expm)
library(matlib)
library(permutations)
library(Rcpp)
sourceCpp("goodness-of-fit-tests\\minor_based_tests_c_optimization.cpp")

#############################
# Generate full minors of M #
#############################
# Calculate all minors of the matrix M consisting of the second and third moments of a 
# distribution to test if a distribution follows a linear SEM corresponding to the full
# graph.

# Calculate all minors of the matrix M for p nodes in the format needed for the not optimized code
calculate_all_minors <- function(p) {
    col_names <- apply(combinations(p, 2, 1:p, repeats.allowed=T), 1, function(row) paste(row, collapse = ""))
    M <- matrix(, nrow = 1 + p, ncol = length(col_names))
    rownames(M) <- c("", 1:p)
    colnames(M) <- col_names
    for (i in 1:nrow(M)) {
        for (j in 1:ncol(M)) {
            if (i == 1) {
                M[i, j] <- paste("s", colnames(M)[[j]], sep = "")
            } else {
                M[i, j] <- paste("t", rownames(M)[[i]], colnames(M)[[j]], sep = "")
            }
        }
    }

    col_subsets <- apply(combinations(ncol(M), p+1, 1:ncol(M), repeats.allowed=F), 1, unlist, simplify = FALSE)
    row_subsets <- apply(combinations(nrow(M), p+1, 1:nrow(M), repeats.allowed=F), 1, unlist, simplify = FALSE)
    all_minors <- list()

    for (row_subset in row_subsets) {
        for (col_subset in col_subsets) {
            all_minors[[length(all_minors) + 1]] <- calculate_determinant(M[row_subset, col_subset])
        }
    }

    return(all_minors)
}

calculate_polynomials_DAG_2_nodes <- function(cause) {
    p <- 2
    col_names <- if (cause == 1) c("11", "12") else c("12", "22")
    M <- matrix(, nrow = 1 + p, ncol = length(col_names))
    rownames(M) <- c("", 1:p)
    colnames(M) <- col_names
    for (i in 1:nrow(M)) {
        for (j in 1:ncol(M)) {
            if (i == 1) {
                M[i, j] <- paste("s", colnames(M)[[j]], sep = "")
            } else {
                M[i, j] <- paste("t", rownames(M)[[i]], colnames(M)[[j]], sep = "")
            }
        }
    }

    col_subsets <- apply(combinations(ncol(M), 2, 1:ncol(M), repeats.allowed=F), 1, unlist, simplify = FALSE)
    row_subsets <- apply(combinations(nrow(M), 2, 1:nrow(M), repeats.allowed=F), 1, unlist, simplify = FALSE)
    all_minors <- list()

    for (row_subset in row_subsets) {
        for (col_subset in col_subsets) {
            all_minors[[length(all_minors) + 1]] <- calculate_determinant(M[row_subset, col_subset])
        }
    }

    return(all_minors)
}

calculate_determinant <- function(A) {
    det <- list()
    n <- nrow(A)
    perms <- allperms(n)
    for (i_perm in seq_len(length(perms))) {
        perm <- perms[i_perm]
        monomial <- c(rep(1, n), sgn(perm))
        vars_in_monomial <- c()
        for (i in seq_len(n))
            vars_in_monomial <- append(vars_in_monomial, A[[i, as.function(perm)(i)]])
        names(monomial) <- c(vars_in_monomial, "coef")
        det[[length(det) + 1]] <- monomial
    }
    return(det)
}

polynomial_sym_rank_3 <- list(
    c(coef=1, t111=1, t222=1, t333=1, t123=1), 
    c(coef=-1, t222=1, t333=1, t112=1, t113=1), 
    c(coef=-1, t333=1, t111=1, t122=1, t223=1), 
    c(coef=-1, t111=1, t222=1, t133=1, t233=1), 
    c(coef=-1, t123=1, t111=1, t223=1, t233=1), 
    c(coef=-1, t123=1, t222=1, t133=1, t113=1), 
    c(coef=-1, t123=1, t333=1, t112=1, t122=1), 
    c(coef=1, t111=1, t122=1, t233=2), 
    c(coef=1, t111=1, t133=1, t223=2), 
    c(coef=1, t222=1, t112=1, t133=2), 
    c(coef=1, t222=1, t233=1, t113=2), 
    c(coef=1, t333=1, t223=1, t112=2), 
    c(coef=1, t333=1, t113=1, t122=2), 
    c(coef=-1, t123=4), 
    c(coef=-2, t123=2, t122=1, t133=1), 
    c(coef=-2, t123=2, t233=1, t112=1), 
    c(coef=-2, t123=2, t113=1, t223=1), 
    c(coef=-3, t123=1, t112=1, t223=1, t133=1), 
    c(coef=-3, t123=1, t113=1, t122=1, t233=1), 
    c(coef=-1, t122=2, t133=2), 
    c(coef=-1, t233=2, t112=2), 
    c(coef=-1, t113=2, t223=2), 
    c(coef=1, t233=1, t112=1, t113=1, t223=1), 
    c(coef=1, t113=1, t223=1, t122=1, t133=1), 
    c(coef=1, t122=1, t133=1, t233=1, t112=1)
)

# Translate from R format to format needed for C++ optimized tests

monomial_to_c_format <- function(monomial) {
  coef <- monomial[["coef"]]
  second_moms <- list()
  third_moms <- list()
  for (var in names(monomial)[names(monomial) != "coef"]) {
    type = substring(var, 1, 1)
    var_deg = monomial[[var]]
    if (type=="s") {
      i = as.numeric(substring(var, 2, 2))
      j = as.numeric(substring(var, 3, 3))
      for (d in 1:var_deg){
        second_moms <- append(second_moms, list(c(i, j)))
      }
    }
    else if (type=="t") {
      i = as.numeric(substring(var, 2, 2))
      j = as.numeric(substring(var, 3, 3))
      k = as.numeric(substring(var, 4, 4))
      for (d in 1:var_deg){
        third_moms <- append(third_moms, list(c(i, j, k)))
      }
    }
  }
  second_moms <- if(length(second_moms) > 0) do.call("rbind", second_moms) else matrix(nrow=0, ncol=0)
  third_moms <- if(length(third_moms) > 0) do.call("rbind", third_moms) else matrix(nrow=0, ncol=0)
  return(list("coef"=coef, "second_moms"=second_moms, "third_moms"=third_moms))
}

polynomial_to_c_format <- function(poly) {
  return(lapply(poly, monomial_to_c_format))
}

###########
# helpers #
###########

permutations <- function(n){
  return(do.call(rbind,combinat::permn(seq(n))))
}

findn <- function(N,D){
  rem <- N%%D
  if (rem==0){
    return(N)
  } else {
    return(N-rem)
  }
}

random_combs <- function(n,k,nr){
  v = seq(n)
  res = matrix(0, nr, k)
  for (i in 1:nr){
    res[i,] = sample(v,k)
  }
  # replacement is possible but very unlikely
  return(res)
}

degree <- function(poly){
    return(max(sapply(poly, function(x){sum(x[1:(length(x)-1)])})))
}

degree_poly_c_format <- function(poly) {
    degrees_monomials = sapply(poly, function(mon) nrow(mon[["second_moms"]]) + nrow(mon[["third_moms"]]))
    return(max(degrees_monomials))
}

create_estimator <- function(poly){

  #Degree of polynomial
  deg = degree(poly)

  # Checks
  if(!is.list(poly)){
    stop("Input polynomial is not in required list format.")
  }
  if(!(deg%%1==0)){
    stop("Degree of polynomial is not an integer.")
  }
  
  # Return estimator function
  return(function(L){
    
    if(!(is.list(L))){
      stop("Input is not a list.")
    }
    if(deg>length(L)){
      stop("Total degree of estimated polynomial has to be smaller or equal than the length of L.")
    }

    poly_est = 0
    for (m in 1:length(poly)){
        monomial = poly[[m]]
        variables = names(monomial)[names(monomial) != "coef"]
        monomial_est = monomial["coef"]
        m_deg = 1
        for (var in variables){
            var_deg = monomial[var]
            type = substring(var, 1, 1)
            if (type=="s"){
                i = as.numeric(substring(var, 2, 2))
                j = as.numeric(substring(var, 3, 3))
                for (d in 1:var_deg){
                  monomial_est = monomial_est * L[[m_deg]][i] * L[[m_deg]][j]
                  m_deg = m_deg + 1
                }
            }
            else if (type=="t"){
                i = as.numeric(substring(var, 2, 2))
                j = as.numeric(substring(var, 3, 3))
                k = as.numeric(substring(var, 4, 4))
                for (d in 1:var_deg){
                  monomial_est = monomial_est * L[[m_deg]][i] * L[[m_deg]][j] * L[[m_deg]][k]
                  m_deg = m_deg + 1
                }
            }
        }
        poly_est = poly_est + monomial_est
    }
    return(as.numeric(poly_est))
  })
}

###########################
## matrix H of estimates ##
###########################

h <- function(L, estimators){
    p = length(estimators)
    h = rep(0,p)
    perm = permutations(length(L))
    for(per in 1:nrow(perm)){
      h = h + sapply(estimators, function(f) f(L[perm[per,]]))
    }
    return(h/nrow(perm))
}

H <- function(X, indices, estimators){
  H = matrix(0,nrow=nrow(indices),ncol=length(estimators))
  for(i in 1:nrow(indices)){
    L = list()
    for (k in 1:ncol(indices)){
      L = c(L, list(X[indices[i,k],]))
    } 
    H[i,] = h(L, estimators)
  }
  return(H)
}

compute_S <- function(n, i, len){
  K = floor((n-1)/len)
  v = rep(0, (n-1))
  for (j in 1:n){
    if (j < i){
      v[j] = j
    } else if (j==i){
      next
    } else {
      v[j-1] = j
    }
  }
  res = matrix(v[1:(K*len)], nrow=K, ncol=len, byrow=TRUE)
  return(res)
}

g <- function(X, i, order_kernel, estimators){
  
  n = nrow(X)
  res = rep(0, length(estimators))
  K = floor((n-1)/(order_kernel-1))
  S = compute_S(n,i,(order_kernel-1))
  
  for (k in 1:K){
    L = list(X[i,])
    for (j in 1:(order_kernel-1)){
      L = c(L, list(X[S[k,j],]))
    }
    res = res + h(L, estimators)
  }
  return(res)
}

G <- function(X, order_kernel, estimators){
  n = nrow(X)
  res = matrix(0, nrow=n, ncol=length(estimators))
  for (i in 1:n){
    res[i,] = g(X, i, order_kernel, estimators)
  }
  return(res)
}

###########
## indep ##
###########

test_indep <- function(X, polynomials, E=1000){

    order_kernel = max(sapply(polynomials, degree))
    N = findn(nrow(X),order_kernel)
    indices = matrix(1:N, ncol=order_kernel)

    # Create estimator functions from polynomials
    estimators <- list()
    for (j in 1:length(polynomials)){
        estimators[[j]] <- create_estimator(polynomials[[j]])
    }

    H = H(X, indices, estimators)
    n = dim(H)[1]
    p = dim(H)[2]  # nr of equality constraints
  
    # Mean and centering
    H_mean = Rfast::colmeans(H)
    H_centered = Rfast::transpose(Rfast::transpose(H) - H_mean) # Centering: H_i = (H_i - H_mean)
    
    # Diagonal of the sample covariance of H
    H_cov = Rfast::colsums(H_centered**2) / n
    
    # Vector for studentizing
    studentizer = H_cov**(-1/2)
    
    # Test statistic
    marginal_stats = studentizer * sqrt(n) * H_mean
    test_stat =  max(abs(marginal_stats))
    
    # Bootstrapping 
    W = matrix(0, E, p)
    for (i in 1:E){
        epsilons = rnorm(n,0,1)
        for (j in 1:p){
            W[i,j] = 1/sqrt(n) * sum(H_centered[,j] * epsilons)
        }
    }
    W = abs(W)
    W_studentized = Rfast::transpose(Rfast::transpose(W) * studentizer)
    results = Rfast::rowMaxs(W_studentized, value = TRUE)
    
    # pval
    pval = (1 + sum(results >= test_stat)) / (1+E)
    
    return(list("PVAL"=pval, "TSTAT"=test_stat, "Results"=results, "H"=H))
}

############
## U-stat ##
############

test_U_stat <- function(X, polynomials, E=1000, n1=min(nrow(X),500), N=2*nrow(X)){
  n = nrow(X)
  p = length(polynomials)
  order_kernel = max(sapply(polynomials, degree_poly_c_format))
  N = min(0.7*choose(n,order_kernel), N)

  # Determine N_hat by Bernoulli sampling
  N_hat = stats::rbinom(1, choose(n,order_kernel), (N / choose(n,order_kernel)))

  # Choose randomly N_hat unique subsets with cardinality r of {1,...,n}
  indices = matrix(unlist(random_combs_c(n,order_kernel,N_hat)[[1]]), ncol = order_kernel, byrow = TRUE)

  # Compute matrix H
  H = H_c(X, indices, polynomials)
  H_mean = colMeans(H)
  H_centered = t(t(H) - H_mean)

  # Compute matrix G on subset of samples
  X_G = X[sample(n, size=min(n1,n), replace=FALSE),]
  G = G_c(X_G, order_kernel, polynomials)
  G_mean = colMeans(G)
  G_centered = t(t(G) - G_mean)

  # Diagonal of the approximate variance of H
  cov_H_diag = colSums(H_centered**2) / N_hat
  cov_G_diag = colSums(G_centered**2) / n1
  cov_diag = order_kernel**2 * cov_G_diag + (n/N) * cov_H_diag

  # Vector for studentizing
  studentizer = cov_diag**(-1/2)

  # Test statistic
  marginal_stats = abs(sqrt(n) * H_mean)
  test_stat =  max(studentizer * marginal_stats)

  # Bootstrap
  U_B = matrix(0, E, p)
  for (i in 1:E){
    epsilons = rnorm(N_hat,0,1)
    for (j in 1:p){
        U_B[i,j] = 1/sqrt(N_hat) * sum(H_centered[,j] * epsilons)
    }
  }
  U_A = matrix(0, E, p)
  for (i in 1:E){
    epsilons = rnorm(n1,0,1)
    for (j in 1:p){
        U_A[i,j] = 1/sqrt(n1) * sum(G_centered[,j] * epsilons)
    }
  }
  U = abs(order_kernel * U_A + sqrt(n/N) * U_B)
  U_studentized = t(t(U) * studentizer)
  results = matrixStats::rowMaxs(U_studentized)

  # pval
  pval = (1 + sum(results >= test_stat)) / (1+E)

  return(list("PVAL"=pval, "TSTAT"=test_stat))
}

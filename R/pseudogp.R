
#' Fit probabilistic pseudotime model
#'
#' @export
fitPseudotime <- function(X, smoothing_alpha = 10, smoothing_beta = 3,
                          pseudotime_mean = 0.5, pseudotime_var = 1,
                          chains = 1, iter = 1000, ...) {
  ## find number of representations
  if(is.matrix(X)) X <- list(X)
  if(!is.list(X)) stop("X must either be matrix (for single representation) or list of matrices")
  if(!all(sapply(X, class) == "matrix")) stop("X must either be matrix (for single representation) or list of matrices")
  dims <- sapply(X, dim)
  if(!all(dims[1,] == dims[1,1]) | !all(dims[2,] == dims[2,1])) stop("All representations must be of same dimensionality")

  ncells <- dim(X[[1]]) [1] # number of cells
  ndim <- dim(X[[1]])[2] # number of dimensions (currently just 2 supported)
  stopifnot(ndim == 2)
  Ns <- length(X) # number of representations

  message(paste("Creating pseudotime model with", ncells, "cells and ", Ns, "representation(s)"))

  ## sanity check (e.g. all representations centred)
  message("Standardizing input data")
  X <- lapply(X, standardize)

  ## prepare data
  #alpha_beta <- meanvar_to_alphabeta(smoothing_mean, smoothing_var)
  #message(paste("Using hyperparameters", alpha_beta))

  # ***always*** order arrays by (representation, latent dimension, cell)
  dx <- array(dim = c(Ns, ndim, ncells))
  for(i in 1:Ns) dx[i,,] <- t(X[[i]])

  data <- list(Ns = Ns, P = ndim, N = ncells,
               X = dx,
               gamma_alpha = smoothing_alpha, gamma_beta = smoothing_beta,
               pseudotime_mean = pseudotime_mean, pseudotime_var = pseudotime_var)

  stanfile <- system.file("pseudogp.stan", package = "pseudogp")

  if(!require(rstan)) stop("Stan required for inference")
  fit <- stan(file = stanfile, data = data,
                    iter = iter, chains = chains, ...)

  return( fit )
}

standardize <- function(X) {
  X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))
}

#' Convert the mean and variance of a gamma distribution
#' to the alpha-beta parametrization
meanvar_to_alphabeta <- function(mean, var) {
  alpha <- mean^2 / var
  beta <- mean / var
  return(c(alpha, beta))
}

#' Calculate a double exponential covariance matrix
#' with a diagonal noise component
cov_matrix <- function(t1, t2, lambda, sigma = NULL) {
  n1 <- length(t1)
  n2 <- length(t2)
  C <- matrix(NA, nrow = n1, ncol = n2)
  for(i in 1:n1) {
    for(j in 1:n2) {
      C[i, j] <- exp(-lambda * (t1[i] - t2[j])^2)
    }
  }
  if(!is.null(sigma)) {
    stopifnot(n1 == n2)
    C <- C + sigma * diag(n1)
  }
  return ( C )
}

#' Compute the posterior mean curve
#'
#' This computes the posterior mean curve. Note that the "mean" aspect here
#' doesn't relate to pseudotimes - this can take a pseudotime sample from the
#' posterior pseudotime distribution and will calculate the mean curve
#' in the data space.
#'
#' @param X The data representation
#' @param t The pseudotime sample
#' @param l Lambda
#' @param s Signa
#' @param nnt Number of new time points
posterior_mean_curve <- function(X, t, l, s, nnt = 80) {
  nt <- runif(nnt)
  K_y <- lapply(1:2, function(i) cov_matrix(t, t, as.numeric(l[i]), as.numeric(s[i])))
  K_star <- lapply(1:2, function(i) cov_matrix(t, nt, as.numeric(l[i])))
  K_dstar <- lapply(1:2, function(i) cov_matrix(nt, nt, as.numeric(l[i])))

  mu_star <- lapply(1:2, function(i) {
    t(K_star[[i]]) %*% solve(K_y[[i]]) %*% X[,i]
  })

  mus <- do.call(cbind, mu_star)
  return(list( mu = mus, t = nt ))
}

#' Laplacian eigenmaps representation of monocle data
"monocle_le"

#' PCA representation of monocle data
"monocle_pca"

#' tSNE representation of monocle data
"monocle_tsne"

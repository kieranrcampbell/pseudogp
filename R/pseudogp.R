
#' Fit probabilistic pseudotime model
#'
#' @export
fitPseudotime <- function(X, smoothing_mean = 6, smoothing_var = 1,
                          pseudotime_mean = 0.5, pseudotime_var = 1,
                          chains = 1, iter = 1000, ...) {
  ## find number of representations
  if(is.matrix(X)) X <- list(X)
  if(!is.list(X)) stop("X must either be matrix (for single representation) or list of matrices")
  if(!all(sapply(X, class) == "matrix")) stop("X must either be matrix (for single representation) or list of matrices")
  dims <- sapply(X, dim)
  if(!all(dims[1,] == dims[1,1]) | !all(dims[2,] == dims[2,2])) stop("All representations must be of same dimensionality")

  ncells <- dim(X[[1]]) [1] # number of cells
  ndim <- dim(X[[1]])[2] # number of dimensions (currently just 2 supported)
  stopifnot(ndim == 2)
  Ns <- length(X) # number of representations

  message(paste("Creating pseudotime model with", ncells, "cells"))

  ## sanity check (e.g. all representations centred)
  message("Standardizing input data")
  X <- lapply(X, standardize)

  ## prepare data
  alpha_beta <- meanvar_to_alphabeta(smoothing_mean, smoothing_var)

  # ***always*** order arrays by (representation, latent dimension, cell)
  dx <- array(dim = c(Ns, ndim, ncells))
  for(i in 1:Ns) dx[i,,] <- X[[i]]

  data <- list(Ns = Ns, P = ndim, N = ncells,
               X = dx,
               gamma_alpha = alpha_beta[1], gamma_beta = alpha_beta[2],
               pseudotime_mean = pseudotime_mean, pseudotime_var = pseudotime_var)

  stanfile <- system.file("pseudogp.stan", package = "pseudogp")

  fit <- rstan::fit(file = stanfile, data = data,
                    iter = iter, chains = chains, ...)
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




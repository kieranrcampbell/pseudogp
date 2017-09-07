
#' Posterior boxplot ordered by pseudotime
#'
#' @param fit The fit object returned by \code{fitPseudotime}
#' @param inner The inner interval for box 'edges' (default 0.75)
#' @param outer The outer interval for boxplot whiskers (default 0.95)
#'
#' @return A \code{ggplot2} boxplot of posterior pseudotime samples ordered by
#' median pseudotime.
#'
#' @export
posteriorBoxplot <- function(fit, inner = 0.75, outer = 0.95) {
  # we can essentially infer the number of chains & representations given the dimension
  # of the initialisation
  chains <- length(fit@inits)
  Ns <- dim(fit@inits[[1]]$lambda)[1]
  P <- dim(fit@inits[[1]]$lambda)[2]
  if(chains > 1) warning("Boxplots currently only supported for 1 chain - will permute samples")

  pst <- rstan::extract(fit, pars = "t", permute = TRUE)$t
  tmcmc <- coda::mcmc(pst)
  hpdinner <- coda::HPDinterval(tmcmc, inner)
  hpdouter <- coda::HPDinterval(tmcmc, outer)
  p <- data.frame(cbind(hpdinner, hpdouter))
  edge_names <- c(paste0(c("lower","upper"), inner), paste0(c("lower","upper"), outer))
  names(p) <- edge_names
  p$Median <- matrixStats::colMedians(pst)
  p$Cell <- as.factor(rank(p$Median))

  ggplot(p) + geom_boxplot(aes_string(x = "Cell", middle = "Median",
                                      lower = edge_names[1], upper = edge_names[2],
                                      ymin = edge_names[3], ymax = edge_names[4]),
                           stat = "identity", fill = "darkred", alpha = 0.5) +
    theme_bw() +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
    xlab("Cell") + ylab("Pseudotime")
}

#' Posterior curve plot.
#'
#' Flexible methods for plotting posterior mean
#' curves or samples of posterior mean curves.
#'
#' @param X The representation(s) passed to \code{fitPseudotime} (either a matrix or
#' list of matrices)
#' @param fit The \code{stanfit} object returned by \code{fitPseudotime}
#' @param posterior_mean Logical. If TRUE (default) then the posterior mean curve is
#' plotting at the MAP estimates of all inferred parameters. If FALSE, then \code{nsamples}
#' will be randomly drawn from the posterior of the parameters and a curve plotted for each.
#' @param nsamples The number of new points at which to calculate the posterior mean curves.
#' @param nnt Number of new pseudo time points at which to plot the posterior mean curves.
#' @param point_colour The colour of the points (cells) to draw
#' @param curve_colour The colour of the curves to draw
#' @param point_alpha The alpha (opacity) of the points
#' @param curve_alpha The alpha (opacity) of the curves. Note that this is just a suggested value
#' and the function will choose an appropriate value depending on the number of samples to plot.
#' This is chosen as \eqn{(1 - \alpha) * exp(1 - nsamples) + \alpha}.
#' @param grid_nrow If more than one representation is present then they're plotted in a grid.
#' By default \code{cowplot} will choose the number of rows, but this overrides.
#' @param grid_nrow If more than one representation is present then they're plotted in a grid.
#' By default \code{cowplot} will choose the number of columns, but this overrides.
#' @param standardize_ranges Logical. If plotting multiple representations it can be useful to have
#' x and y lims that don't depend on the fit (so plots align correctly). If this is set to FALSE,
#' \code{ggplot2} calculates the x and y limits. If this is set to TRUE, the x and y limits are set
#' to the minimum and maximum of the X values plus or minus 6\% of the range between them.
#'
#' @export
posteriorCurvePlot <- function(X, fit, posterior_mean = TRUE,
                               nsamples = 50, nnt = 80,
                               point_colour = "darkred",
                               curve_colour = "black", point_alpha = 1,
                               curve_alpha = 0.5,
                               grid_nrow = NULL, grid_ncol = NULL,
                               use_cowplot = TRUE,
                               standardize_ranges = FALSE, ...) {
  if(is.matrix(X)) X <- list(X)
  Ns <- length(X) ## number of representations
  chains <- length(fit@inits)
  message(paste("Plotting traces for", Ns,"representation(s) and", chains, "chain(s)"))

  plots <- vector("list", Ns)
  # this is of dim trace-chain-cell
  pst <- rstan::extract(fit, pars = "t", permute = FALSE)
  lambda <- rstan::extract(fit, pars = "lambda", permute = FALSE)
  sigma <- rstan::extract(fit, pars = "sigma", permute = FALSE)
  for(i in 1:Ns) {
    l <- lambda[,,(2*i - 1):(2*i),drop=FALSE]
    s <- sigma[,,(2*i - 1):(2*i),drop=FALSE]
    plt <- makeEnvelopePlot(pst, l, s, X[[i]], chains, posterior_mean, nsamples, nnt, point_colour,
                            curve_colour, point_alpha, curve_alpha, use_cowplot, standardize_ranges)
    plots[[i]] <- plt
  }
  gplt <- cowplot::plot_grid(plotlist = plots, labels = names(X), nrow = grid_nrow, ncol = grid_ncol)
  return( gplt )
}

#' @importFrom MCMCglmm posterior.mode
makeEnvelopePlot <- function(pst, l, s, x, chains, posterior_mean, ncurves, nnt,
                             point_colour = "darkred", curve_colour = "black",
                             point_alpha = 1, curve_alpha = 0.5,
                             use_cowplot = TRUE, standardize_ranges = FALSE) {
  n_posterior_samples <- dim(pst)[1]
  curve_samples <- sample(n_posterior_samples, ncurves)
  pmcs <- lapply(1:chains, function(chain) {
    if(posterior_mean) {
      tmap <- MCMCglmm::posterior.mode(coda::mcmc(pst[,chain,]))
      lmap <- MCMCglmm::posterior.mode(coda::mcmc(l[,chain,]))
      smap <- MCMCglmm::posterior.mode(coda::mcmc(s[,chain,]))
      return( list( posterior_mean_curve(x, tmap, lmap, smap, nnt) ) )
    } else {
      lapply(curve_samples, function(i) {
        t <- pst[i,chain,]
        lambda <- l[i,chain,]
        sigma <- s[i,chain,]
        posterior_mean_curve(x, t, l, s, nnt)
      })
    }
  })

  x <- as.data.frame(x)
  names(x) <- c("x1", "x2")
  plt <- ggplot()
  plt <- plt + geom_point(data = data.frame(x), aes(x = x1, y = x2), shape = 21,
                          fill = point_colour, colour = 'white', size = 3, alpha = point_alpha) +
    xlab("Component 1") + ylab("Component 2")
  if(standardize_ranges) { # hard set x and y lims
    mins <- apply(x, 2, min)
    maxs <- apply(x, 2, max)
    ranges <- maxs - mins
    pct_range_6 <- 0.06 * ranges
    lower_lims <- mins - pct_range_6
    upper_lims <- maxs + pct_range_6
    plt <- plt + xlim(c(lower_lims[1], upper_lims[1])) + ylim(c(lower_lims[2], upper_lims[2]))
  }

  if(use_cowplot) {
    plt <- plt + cowplot::theme_cowplot()
  } else {
    plt <- plt + theme_bw()
  }

#   ncolor <- min(chains, 9)
#   if(ncolor < 3) ncolor <- 3
#   getPalette <- colorRampPalette(brewer.pal(ncolor, "Set1"))
#   colorset <-getPalette(chains)

  for(chain in 1:chains) {
    pmc <- pmcs[[chain]]
    ncurves <- length(pmc)
    mus <- lapply(pmc, function(p) p$mu)
    M <- data.frame(do.call("rbind", mus))
    names(M) <- c("M1", "M2")
    M$curve <- rep(1:ncurves, each = nrow(mus[[1]]))
    M$nt <- unlist(lapply(pmc, function(x) x$t))
    M <- dplyr::arrange(M, curve, nt)

    calculated_alpha <- (1 - curve_alpha) * exp(1) * exp(-ncurves) + curve_alpha
    if(posterior_mean) calculated_alpha <- 1

    for(i in 1:ncurves) {
      plt <- plt + geom_path(data = dplyr::filter(M, curve == i), aes(x = M1, y = M2),
                             size = 2, alpha = calculated_alpha, color = curve_colour)
    }
  }
  return( plt )
}

#' Plot MCMC diagnostics.
#'
#' Plot basic MCMC diagnostics (traceplot and autocorrelation) of the log-posterior probability
#' for a \code{stanfit} object.
#'
#' Further assessment of convergence can be done using \code{rstan} functions.
#'
#' @param fit A \code{stanfit} object
#' @param arrange How to arrange the plots. If "vertical", traceplot and autocorrelation are
#' arranged in one column, while if "horizontal" traceplot and autocorrelation are arranged
#' in one row.
#' @export
#' @importFrom cowplot plot_grid
#'
#' @return A \code{ggplot2} object
#'
plotDiagnostic <- function(fit, arrange = c("vertical", "horizontal")) {
  stopifnot(is(fit, "stanfit"))
  arrange <- match.arg(arrange)
  nrow <- switch(arrange,
                 vertical = 2,
                 horizontal = 1)
  plt <- cowplot::plot_grid(stan_trace(fit, "lp__"), stan_ac(fit, "lp__"), nrow = nrow)
  return(plt)
}

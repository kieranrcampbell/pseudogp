
#' Posterior boxplot ordered by pseudotime
#'
#' @param fit The fit object returned by \code{fitPseudotime}
#' @param inner The inner interval for box 'edges' (default 0.75)
#' @param outer The outer interval for boxplot whiskers (default 0.95)
#'
#' @export
posteriorBoxplot <- function(fit, inner = 0.75, outer = 0.95) {
  # we can essentially infer the number of chains & representations given the dimension
  # of the initialisation
  chains <- length(fit@inits)
  Ns <- dim(fit@inits[[1]]$lambda)[1]
  P <- dim(fit@inits[[1]]$lambda)[2]
  if(chains > 1) warning("Boxplots currently only supported for 1 chain - will permute samples")

  pst <- extract(fit, pars = "t", permute = TRUE)$t
  tmcmc <- mcmc(pst)
  hpdinner <- HPDinterval(tmcmc, inner)
  hpdouter <- HPDinterval(tmcmc, outer)
  p <- data.frame(cbind(hpdinner, hpdouter))
  edge_names <- c(paste0(c("lower","upper"), inner), paste0(c("lower","upper"), outer))
  names(p) <- edge_names
  p$Median <- colMedians(pst)
  p$Cell <- as.factor(rank(p$Median))

  ggplot(p) + geom_boxplot(aes_string(x = "Cell", middle = "Median",
                                      lower = edge_names[1], upper = edge_names[2],
                                      ymin = edge_names[3], ymax = edge_names[4]),
                           stat = "identity", fill = "darkred", alpha = 0.5) +
    theme_bw() +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
    xlab("Cell") + ylab("Pseudotime")
}

#' Posterior curve plot (workhorse of the analysis)
#'
#' @export
posteriorCurvePlot <- function(X, fit, nsamples = 50, nnt = 80, ...) {
  Ns <- length(X) ## number of representations
  chains <- length(fit@inits)
  message(paste("Plotting traces for", Ns,"representations and", chains, "chains"))

  plots <- vector("list", Ns)
  # this is of dim trace-chain-cell
  pst <- extract(fit, pars = "t", permute = FALSE)
  lambda <- extract(fit, pars = "lambda", permute = FALSE)
  sigma <- extract(fit, pars = "sigma", permute = FALSE)
  for(i in 1:Ns) {
    l <- lambda[,,(2*i - 1):(2*i)]
    s <- sigma[,,(2*i - 1):(2*i)]
    plt <- makeEnvelopePlot(pst, l, s, X[[i]], chains, nsamples, nnt)
    plots[[i]] <- plt
  }
  gplt <- cowplot::plot_grid(plotlist = plots, labels = names(X))
  return( gplt )
}

makeEnvelopePlot <- function(pst, l, s, x, chains, ncurves, nnt) {
  n_posterior_samples <- dim(pst)[1]
  curve_samples <- sample(n_posterior_samples, ncurves)
  pmcs <- lapply(1:chains, function(chain) {
    lapply(curve_samples, function(i) {
      t <- pst[i,chain,]
      lambda <- l[i,chain,]
      sigma <- s[i,chain,]
      posterior_mean_curve(x, t, l, s, nnt)
    })
  })

  x <- as.data.frame(x)
  names(x) <- c("x1", "x2")
  plt <- ggplot()
  plt <- plt + geom_point(data = data.frame(x), aes(x = x1, y = x2), shape = 21,
                          fill = 'black', colour = 'white', size = 3, alpha = 0.5) +
    xlab("") + ylab("")

  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  colorset <-getPalette(chains + 1)

  for(chain in 1:chains) {
    pmc <- pmcs[[chain]]
    ncurves <- length(pmc)
    mus <- lapply(pmc, function(p) p$mu)
    M <- data.frame(do.call("rbind", mus))
    names(M) <- c("M1", "M2")
    M$curve <- rep(1:ncurves, each = nrow(mus[[1]]))
    M$nt <- unlist(lapply(pmc, function(x) x$t))
    M <- dplyr::arrange(M, curve, nt)

    for(i in 1:ncurves) {
      plt <- plt + geom_path(data = dplyr::filter(M, curve == i), aes(x = M1, y = M2),
                             size = 2, alpha = .2, color = colorset[chain])
    }
  }
  return( plt )
}

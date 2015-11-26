
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
posteriorCurvePlot <- function(X, fit, ...) {

}

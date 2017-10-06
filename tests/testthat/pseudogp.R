context("Pseudogp")

test_that("Pseudogp works", {
  x <- runif(100, -1, 1)
  y <- x^2 + runif(100, -.1, .1)
  space <- standardize(cbind(x, y))

  fit <- fitPseudotime(space)

  pdf("/dev/null")
  posteriorCurvePlot(space, fit)
  posteriorBoxplot(fit)
  dev.off()
})

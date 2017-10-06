context("Pseudogp")

test_that("Pseudogp works", {
  x <- runif(100, -1, 1)
  y <- x^2 + runif(100, -.1, .1)
  space <- standardize(cbind(x, y))

  # do not expect error
  expect_error({
    fit <- fitPseudotime(space)
  }, NA)

  expect_is(fit, "stanfit")

  g1 <- posteriorCurvePlot(space, fit)
  g2 <- posteriorBoxplot(fit)

  expect_is(g1, "ggplot")
  expect_is(g2, "ggplot")

  # expect no errors
  expect_error({
    pdf("/dev/null")
    print(g1)
    print(g2)
    dev.off()
  }, NA)
})

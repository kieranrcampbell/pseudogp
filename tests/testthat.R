library(testthat)

# see https://github.com/hadley/testthat/issues/86
Sys.setenv("R_TESTS" = "")

test_check("pseudogp")

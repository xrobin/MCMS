library(testthat)
library(xavamess)

Sys.unsetenv("R_TESTS")
test_check("xavamess")

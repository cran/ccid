context("The change-points in the mean structure for a multivariate
        time series of multiplicative structure via minimising an information
        criterion related to the L2 distance")

# In this file we test the ic_L2 function

test_that("The correct change-points are given when the function is used correctly", {
  set.seed(1)
  num.nodes <- 40 # number of nodes
  etaA.1    <- 0.95
  etaA.2    <- 0.05
  pcor1     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.1)
  pcor2     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.2)

  n <- 100
  data1 <- GeneNet::ggm.simulate.data(n, pcor1)
  data2 <- GeneNet::ggm.simulate.data(n, pcor2)

  X <- rbind(data1, data2) ## change-point at 100

  sgn <- sign(stats::cor(X))

  expect_equal(ic_L2(t(hdbinseg::gen.input(t(X), -1, TRUE, diag = FALSE, sgn)))$changepoints, c(100))
})

test_that("Error and warning messages are given correctly", {
  set.seed(100)
  test_data <- matrix(data = rchisq(1000, df = 5), nrow = 500, ncol = 2)

  # A string vector is given as the input for X.
  expect_error(ic_L2(X = "I am not a numeric matrix"))

  # The threshold constant is set equal to zero.
  expect_error(ic_L2(X = test_data, thr_const = 0))

  # The value for lambda is negative.
  expect_error(ic_L2(X = test_data, points = -2))

  # The value for lambda is a positive real number instead of a positive integer.
  expect_warning(ic_L2(X = test_data, points = 3.4))
})

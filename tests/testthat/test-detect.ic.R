context("The change-points, through model optimisation, in the cross-covariance
        structure for a multivariate time series")

# In this file we test the detect.ic function

test_that("The correct result are given when the function is used correctly", {
  set.seed(111)
  num.nodes <- 40 # number of nodes
  etaA.1    <- 0.95
  etaA.2    <- 0.05
  pcor1     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.1)
  pcor2     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.2)

  n <- 100
  data1 <- GeneNet::ggm.simulate.data(n, pcor1)
  data2 <- GeneNet::ggm.simulate.data(n, pcor2)
  X1 <- rbind(data1, data2) ## change-point at 100

  expect_equal(detect.ic(X1, approach = "euclidean")$changepoints, c(100))
  expect_equal(detect.ic(X1, approach = "infinity")$changepoints, c(100))
})

test_that("Error and warning messages are given correctly", {
  set.seed(100)
  test_data <- matrix(data = rchisq(50000, df = 5), nrow = 500, ncol = 100)

  # A string vector is given as the input for X.
  expect_error(detect.ic(X = "I am not a numeric matrix"))

  # The threshold constant is set equal to zero.
  expect_error(detect.ic(X = test_data, approach = "infinity", th_max = -5))
  expect_error(detect.ic(X = test_data, th_sum = -5))

  # The value for lambda is negative.
  expect_error(detect.ic(X = test_data, pointsgen = -2))

  # The value for lambda is a positive real number instead of a positive integer.
  expect_warning(detect.ic(X = test_data, pointsgen = 3.4))
})

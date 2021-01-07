context("The change-points in the mean structure for a multivariate time series
        of multiplicative structure through thresholding and extraction of the
        component time series where the changes occurred via taking the square
        root of the average of squares of contrast function values")

# In this file we test the cpt_ts_L2 function

test_that("The correct results are given when the function is used correctly", {
  set.seed(1)
  num.nodes <- 40 # number of nodes
  etaA.1    <- 0.95
  etaA.2    <- 0.05
  pcor1     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.1)
  pcor2     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.2)

  n <- 100
  data1 <- GeneNet::ggm.simulate.data(n, pcor1)
  data2 <- GeneNet::ggm.simulate.data(n, pcor2)

  X <- rbind(data1, data2, data1, data2) ## change-points at 100, 200, 300
  sgn <- sign(stats::cor(X))

  expect_equal(cpt_ts_L2(t(hdbinseg::gen.input(t(X), -1, TRUE, diag = FALSE, sgn)))$changepoints, c(100, 200, 300))
  expect_true(is.list(cpt_ts_L2(t(hdbinseg::gen.input(t(X), -1, TRUE, diag = FALSE, sgn)))$time_series))
})

test_that("Error and warning messages are given correctly", {
  set.seed(100)
  test_data <- matrix(data = rchisq(63750, df = 5), nrow = 50, ncol = 1275)

  # A string vector is given as the input for X.
  expect_error(cpt_ts_L2(X = "I am not a numeric matrix"))

  # The threshold constant is set equal to zero.
  expect_error(cpt_ts_L2(X = test_data, thr_const_s = 0))

  # The value for lambda is negative.
  expect_error(cpt_ts_L2(X = test_data, points_s = -2))

  # The value for lambda is a positive real number instead of a positive integer.
  expect_warning(cpt_ts_L2(X = test_data, points_s = 3.4))
})

context("The function that preaverages the given data")

# In this file we test the preaverage function

test_that("The function gives reasonable results when used correctly by the user", {
  expect_equal(dim(preaverage(matrix(rchisq(1000, df = 5), 100, 10), scal = 3)), c(34, 10))
  expect_equal(dim(preaverage(matrix(rchisq(1000, df = 5), 100, 10), scal = 5)), c(20, 10))
})

test_that("The function gives error messages whenever the input argument is meaningless", {
  set.seed(100)
  test_data <- matrix(data = rchisq(1000, df = 5), nrow = 100, ncol = 10)

  # A string vector is given as the input for X.
  expect_error(preaverage(X = "I am not a numeric matrix"))

  # The value for the scaling parameter scal should be positive.
  expect_error(preaverage(X = test_data, scal = -5))

  # The value for the scaling parameter scal is a positive real number instead of a positive integer.
  expect_warning(preaverage(X = test_data, scal = 3.4))
})

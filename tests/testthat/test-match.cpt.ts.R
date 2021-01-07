context("We perform a contrast function based approach in order to give
        information on the time series that contain each change-point. In simple
        words, for a given change-point set this function associates each
        change-point with the respective data sequence (or sequences) from which
        it was detected.")

# In this file we test the match.cpt.ts function

test_that("The correct result is given when the function is used correctly", {
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
  expect_true(is.list(match.cpt.ts(t(hdbinseg::gen.input(t(X), -1, TRUE, diag = FALSE, sgn)), cpt = c(100, 200, 300))))
})

test_that("The function still works if the user does not give a vector of change-points", {
  set.seed(1)
  num.nodes <- 40 # number of nodes
  etaA.1    <- 0.95
  etaA.2    <- 0.05
  pcor1     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.1)
  pcor2     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.2)

  n <- 100
  data1 <- GeneNet::ggm.simulate.data(n, pcor1)
  data2 <- GeneNet::ggm.simulate.data(n, pcor2)
  test_data <- rbind(data1, data2, data1, data2) ## change-points at 100, 200, 300

  sgn <- sign(cor(test_data))
  expect_true(is.list(match.cpt.ts(t(hdbinseg::gen.input(t(test_data), -1, TRUE, diag = FALSE, sgn)))))
})

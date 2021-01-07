
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ccid

<!-- badges: start -->

<!-- badges: end -->

The goal of ccid is to implements the Cross-Covariance Isolate Detect
(CCID) methodology for the estimation of the number and location of
multiple change-points in the second-order (cross-covariance or network)
structure of multivariate, possibly high-dimensional time series. The
method is motivated by the detection of change points in functional
connectivity networks for functional magnetic resonance imaging (fMRI),
electroencephalography (EEG), magentoencephalography (MEG) and
electrocorticography (ECoG) data. The main routines in the package have
been extensively tested on fMRI data. For details on the CCID
methodology, please see Anastasiou et al (2020).

## Installation

You can install the released version of ccid from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ccid")
```

## Example

These are two basic examples which show you how to detect changes (if
there are any) in the second-order (cross-covariance or network)
structure of multivariate, possibly high-dimensional time series.

``` r
library(ccid)
## An example of three change-points in the cross-covariance structure
## of a multivariate time series of length 400 and dimensionality equal to 40.
set.seed(111111)
num.nodes <- 40 # number of nodes
etaA.1    <- 0.95
etaA.2    <- 0.05
pcor1     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.1)
pcor2     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.2)

n <- 100
data1 <- GeneNet::ggm.simulate.data(n, pcor1)
data2 <- GeneNet::ggm.simulate.data(n, pcor2)

X1 <- rbind(data1, data2, data1, data2) ## change-points at 100, 200, 300
N1 <- detect.ic(X1, approach = 'euclidean', scales = -1)
N2 <- detect.ic(X1, approach = 'infinity', scales = -1)
N1$changepoints
#> [1] 100 199 300
N2$changepoints
#> [1] 100 199 300
N1$no.of.cpts
#> [1] 3
N2$no.of.cpts
#> [1] 3

## An example of no change-points.
set.seed(11)
A <- matrix(rnorm(20*400), nrow = 400) ## No change-point
M1 <- detect.ic(A, approach = 'euclidean', scales = -1)
M2 <- detect.ic(A, approach = 'infinity', scales = -1)
M1$changepoints
#> [1] NA
M2$changepoints
#> [1] NA
```

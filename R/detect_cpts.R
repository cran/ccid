options(expressions = 5e+05)
#' Preaveraging the multivariate time series
#'
#' This function pre-processes the given data in order to remove serial
#' correlation that might exist in the given data.
#'
#' @export
#' @param X A numerical matrix representing the multivariate time series,
#'   with the columns representing its components.
#' @param scal A positive integer number with default value equal to 3. It is
#'   used to define the way we pre-average the data sequences.
#' @details For a given natural number \code{scal} and data matrix \code{X} of
#' dimensionality \eqn{T \times d}, let us denote by
#' \eqn{Q = \lceil T/scal \rceil}. Then, \code{\link{preaverage}} calculates,
#' for all \eqn{j = 1,2, ..., d},
#' \deqn{\tilde{X}_{q, j} = 1/scal\sum_{t=(q-1) * sc + 1}^{q * sc}X_{t, j},}
#' for \eqn{q=1, 2, ..., Q-1}, while
#' \deqn{\tilde{x}_{Q, j} = (T - (Q-1) * sc)^{-1}\sum_{t = (Q-1) * sc + 1}^{T}X_{t, j}.}
#' @return
#' The ``preaveraged'' matrix \eqn{\tilde{X}} of dimensionality
#' \eqn{Q \times d}, as explained in Details.
#' @author Andreas Anastasiou, \email{anastasiou.andreas@ucy.ac.cy}
#' @references ``Cross-covariance isolate detect: a new change-point method
#'   for estimating dynamic functional connectivity'', Anastasiou et al (2020),
#'   preprint <doi:10.1101/2020.12.20.423696>.
#' @examples
#' A <- matrix(1:32, 8, 4)
#' A
#' A1 <- preaverage(A, scal = 3)
#' A1
preaverage <- function(X, scal = 3) {
    if (!(is.numeric(X))) {
        stop("The input in `X' should be a numeric matrix which has at each
        column the time series that we want to pre-average.")
    }
    if ((scal <= 0)) {
        stop("The scaling constant which is used to define the way we
        preaverage the data should be a positive number.")
    }
    if (abs(scal - round(scal)) > .Machine$double.eps^0.5) {
        warning("The input for `scal' should be a positive integer. If it is
        a positive real number then the integer part of the given number is
        used as the value of `scal'.")
    }
    scal <- as.integer(scal)
    l <- nrow(X)
    d <- ncol(X)
    l1 <- ceiling(l / scal)
    res <- matrix(NA, l1, d)
    if (l1 == 1) {
        res <- apply(X[((l1 - 1) * scal + 1):l, , drop = FALSE], 2, mean)
    } else {
        for (i in 1:(l1 - 1)) {
            res[i, ] <- apply(X[((i - 1) * scal + 1):(i * scal), , drop = FALSE], 2, mean)
        }
        res[l1, ] <- apply(X[((l1 - 1) * scal + 1):l, , drop = FALSE], 2, mean)
    }
    return(res)
}

#' Associating the change-points with the component time series
#'
#' This function performs a contrast function based approach in order to
#' match each change-point and time series. In simple terms, for a given
#' change-point set this function associates each change-point with the
#' respective data sequence (or sequences) from which it was detected.
#'
#' @export
#' @param X A numerical matrix representing the multivariate periodograms.
#'   Each column contains a different periodogram which is the result of
#'   applying the wavelet transformation to the initial multivariate time
#'   series.
#' @param cpt A positive integer vector with the locations of the
#'   change-points. If missing, then our approach with the \eqn{L_2}
#'   aggregation is called internally to extract the change-points in
#'   \code{X}.
#' @param thr_const A positive real number with default value equal to 1. It is
#'   used to define the threshold; see \code{thr_fin}.
#' @param thr_fin With \code{T} the length of the data sequence, this is a
#'   positive real number with default value equal to
#'   \code{thr_const * log(T)}. It is the threshold, which is used in the
#'   detection process.
#' @param scales Negative integers for the wavelet scales used to create the periodograms,
#'   with a small negative integer representing a fine scale. The default value is equal
#'   to -1.
#' @param count Positive integer with default value equal to 5. It can be used
#'   so that the function will return only the \code{count} most important
#'   matches of each change-points with the time series.
#' @return
#' A list with the following components:
#'   \tabular{ll}{
#'    \cr \code{time_series_indicator} \tab A list of matrices. There are as many matrices as
#'    \cr \tab the number of change-points. Each change-point has its own matrix, with
#'    \cr \tab each row of the matrix representing the associated combination of time
#'    \cr \tab series that are associated with the respective change-point.
#'    \cr \code{most_important} \tab A list of matrices. There are as many matrices as
#'    \cr \tab the number of change-points. Each change-point has its own matrix, with
#'    \cr \tab each row of the matrix representing the associated combination of time
#'    \cr \tab series that are associated with the respective change-point. It shows the
#'    \cr \tab \code{count} most important time series combinations for each change-point.
#'  }
#' @author Andreas Anastasiou, \email{anastasiou.andreas@ucy.ac.cy}
#' @references ``Cross-covariance isolate detect: a new change-point method
#'   for estimating dynamic functional connectivity'', Anastasiou et al (2020),
#'   preprint <doi:10.1101/2020.12.20.423696>.
#' @examples
#'   set.seed(1)
#'   num.nodes <- 40 # number of nodes
#'   etaA.1    <- 0.95
#'   etaA.2    <- 0.05
#'   pcor1     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.1)
#'   pcor2     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.2)
#'
#'   n <- 100
#'   data1 <- GeneNet::ggm.simulate.data(n, pcor1)
#'   data2 <- GeneNet::ggm.simulate.data(n, pcor2)
#'   X <- rbind(data1, data2, data1, data2) ## change-points at 100, 200, 300
#'   sgn <- sign(stats::cor(X))
#'   M1 <- match.cpt.ts(t(hdbinseg::gen.input(x = t(X),scales = -1, sq = TRUE,
#'   diag = FALSE, sgn = sgn)))
#'   M1
match.cpt.ts <- function(X, cpt, thr_const = 1,
                         thr_fin = thr_const * sqrt(2 * log(nrow(X))),
                         scales = -1, count = 5) {
  if (missing(cpt)) {
    cpt <- thr_L2(X)
  }
  lsc <- length(scales)
  s_cpt <- sort(cpt)
  lx <- nrow(X)
  lc <- ncol(X)
  lc_i <- lc/lsc
  count <- min(lc_i, count)
  pdim <- (-lsc + sqrt((lsc^2)+8*lc*lsc))/(2*lsc)
  if (round(pdim) != pdim){
    stop("The input in X of the function match.cpt.ts does not seem to be
         the result of a wavelet transformation on the original time series since
         the dimensionality does not match")
  }
  lcpt <- length(cpt)
  if (lcpt == 0){
    return(list(time_series_indicator = "Matching the time series with the change-points is not possible when no change-points are given."))
  }
  else{
    seb_set <- unique(c(1, s_cpt, lx))
    lseb_set <- length(seb_set)
    Res <- list()
    Res$matches <- list()
    Res$most.important <- list()
    for (j in 1:lsc){
    Res$matches[[j]] <- list()
    Res$most.important[[j]] <- list()
    CS <- matrix(NA, lc_i, lcpt)
    ts1 <- matrix(0, lc_i, 2*lcpt)
    ts1[,seq(1, 2*lcpt, 2)] <- seq(1:(lc_i))
    X_temp <- X[,((j-1)*lc_i + 1) : (j*lc_i)]
    for (i in 1:lc_i) {
      CS[i, ] <- cusum_one_multi(X_temp[,i], seb_set[1:(lseb_set - 2)], seb_set[3:(lseb_set)], seb_set[2:(lseb_set - 1)])
    }
    ts1[, seq(2,2*lcpt,2)] <- CS[, 1:lcpt]
    ord <- list()
    for(i in 1:lcpt){
      ord[[i]] <- order(CS[,i], decreasing=T)
      ts1[,(2*i - 1):(2*i)] <- ts1[ord[[i]],(2*i - 1):(2*i)]
    }
    indic2 <- integer()
    for (i in 1:(lcpt)) {
      indic2[i] <- which(ts1[,2*i] <= thr_fin)[1]
    }
    mind <- max(indic2)
    if (mind == 1){
      ts <- matrix(0,0,2*lcpt)
    }
    else{
      ts <- ts1[1:(mind - 1), ,drop = FALSE]
    }
    res <- list()
    A <- matrix(0, pdim, pdim)
    gdata::upperTriangle(A, diag = TRUE, byrow = T) <- 1:lc_i
    B <- list()
    B$matches <- list()
    for (i in 1:lcpt){
      res[[i]] <- list()
      res[[i]] <- ts[which(ts[,2*i] != 0),(2*i - 1)]
      M_temp <- matrix(0,count,2)
      for(k in 1:count){
        M_temp[k,] <- which(A == ts1[k,(2*i -1)], arr.ind = TRUE)
      }
      M_temp <- cbind(M_temp,ts1[1:count,(2*i)])
      B$most.important[[i]] <- M_temp
      if (length(res[[i]]) == 0) {
        B$matches[[i]] <- matrix(0,0,2)
      }
      else{
        B$matches[[i]] <- matrix(0,length(res[[i]]),2)
        for(k in 1:length(res[[i]])){
          B$matches[[i]][k,] <- which(A == res[[i]][k], arr.ind = TRUE)
        }
      }
    }
    Res$matches[[j]] <- B$matches
    Res$most.important[[j]] <- B$most.important
    }
    Res_fin <- list()
    Res_fin$matches <- list()
    Res_fin$most.important <- list()
    for (i in 1:lcpt){
      Res_fin$matches[[i]] <- matrix(0,1,2)
      Res_fin$most.important[[i]] <- matrix(0,1,3)
    for (j in 1:lsc){
      Res_fin$matches[[i]] <- rbind(Res_fin$matches[[i]], Res$matches[[j]][[i]])
      Res_fin$most.important[[i]] <- rbind(Res_fin$most.important[[i]], Res$most.important[[j]][[i]])
    }
      Res_fin$matches[[i]] <- Res_fin$matches[[i]][-1,]
      Res_fin$most.important[[i]] <- Res_fin$most.important[[i]][-1,]
      Res_fin$matches[[i]] <- unique(Res_fin$matches[[i]])
      ord.temp <- order(Res_fin$most.important[[i]][,3], decreasing=T)
      Res_fin$most.important[[i]] <- unique(Res_fin$most.important[[i]][ord.temp,1:2])
      Res_fin$most.important[[i]] <- Res_fin$most.important[[i]][1:count,]
    }
    return(list(time_series_indicator = Res_fin$matches, most_important = Res_fin$most.important))
  }
}

#' Multiple change-point detection in the cross-covariance structure
#' of multivariate high-dimensional time series using a thresholding
#' based procedure and, wherever possible, extraction of the component
#' time series where the changes occurred
#'
#' This function detects multiple change-points in the cross-covariance
#' structure of a multivariate time series using a thresholding based procedure.
#' It also, wherever possible, returns the relevant, transformed
#' time series where each change-point was detected. See Details for a brief
#' explanation.
#'
#' @export
#' @param X A numerical matrix representing the multivariate time series,
#'   with the columns representing its components.
#' @param approach A character string, which defines the metric to be used
#'   in order to detect the change-points. If approach = ``euclidean'', which
#'   is also the default value, then the \eqn{L_2} metric will be followed
#'   for the detection. If approach = ``infinity'', then the \eqn{L_{\infty}}
#'   metric will be used for the detection.
#' @param th_max A positive real number with default value equal to 2.25. It is
#'   used to define the threshold if the \eqn{L_{\infty}} metric is chosen in
#'   \code{approach} .
#' @param th_sum A positive real number with default value equal to 0.65. It is
#'   used to define the threshold if the \eqn{L_2} metric is chosen in
#'   \code{approach}.
#' @param pointsgen A positive integer with default value equal to 10. It
#'   defines the distance between two consecutive end- or start-points of
#'   the right- or left-expanding intervals, respectively; see Details for
#'   more information.
#' @param scales Negative integers for wavelet scales, with a small negative
#'   integer representing  a fine scale. The default value is equal to -1.
#' @param preaverage_gen A logical variable with default value equal to
#'   \code{FALSE}. If \code{FALSE}, then pre-averaging the data is not
#'   required. If \code{TRUE}, then we need to pre-average the data before
#'   proceeding with the detection of the change-points.
#' @param scal_gen A positive integer number with default value equal to 3. It
#'   is used to define the way we pre-average the given data sequence only if
#'   \code{preaverage_gen = TRUE}. See the Details in \code{\link{preaverage}}
#'   for more information on how we pre-average.
#' @param min_dist A positive integer number with default value equal to 1. It
#'   is used in order to provide the minimum distance acceptable between
#'   detected change-points if such restrictions apply.
#' @details The time series \eqn{X_t} is of dimensionality \eqn{p} and we are
#'   looking for changes in the cross-covariance structure between the
#'   different time series components
#'   \eqn{X_{t}^{(1)}, X_{t}^{(2)}, ..., X_{t}^{(p)}}. We first use a
#'   wavelet-based approach for the various given scales in \code{scales}
#'   in order to transform the given time series \eqn{X_t} to a multiplicative
#'   model \eqn{Y_{t}^{(k)} = \sigma^{(k)}_t (Z_t^{(k)})^2; t=1,2,\ldots,T; k = 1,2,\ldots,d,}
#'   where \eqn{Z_t^{(k)}} is a sequence of standard normal random variables,
#'   \eqn{E(Y_t^{(k)}) = \sigma_t^{(k)}}, and \eqn{d} is the new
#'   dimensionality, which depends on the value given in \code{scales}.
#'   The function has been extensively tested on fMRI data, hence, its parameters
#'   have been tuned for this data type. The function might not work well in other
#'   structures, such as time series that are negatively serially correlated.
#' @return
#' A list with the following components:
#' \tabular{ll}{
#' \cr \code{changepoints} \tab The locations of the detected change-points.
#' \cr \code{no.of.cpts} \tab The number of the detected change-points.
#' \cr \code{time_series} \tab A list with two components that indicates which combinations
#' \cr \tab of time series are responsible for each change-point detected. See the outcome
#' \cr \tab values \code{time_series_indicator} and \code{most_important} of the function
#' \cr \tab \code{\link{match.cpt.ts}} for more information.
#'    }
#'   If the minimum distance between the detected change-points is less than
#'   the value given in \code{min_dist}, then only the number and the locations of
#'   the ``pruned'' change-points are returned.
#' @author Andreas Anastasiou, \email{anastasiou.andreas@ucy.ac.cy}
#' @references ``Cross-covariance isolate detect: a new change-point method
#'   for estimating dynamic functional connectivity'', Anastasiou et al (2020),
#'   preprint <doi:10.1101/2020.12.20.423696>.
#' @seealso \code{\link{detect.ic}}.
#' @examples
#'   set.seed(111)
#'   A <- matrix(rnorm(20*400), nrow = 400) ## No change-point
#'   M1 <- detect.th(A, approach = 'euclidean', scales = -1)
#'   M2 <- detect.th(A, approach = 'infinity', scales = -1)
#'   M1
#'   M2
#'
#'   set.seed(111)
#'   num.nodes <- 40 # number of nodes
#'   etaA.1    <- 0.95
#'   etaA.2    <- 0.05
#'   pcor1     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.1)
#'   pcor2     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.2)
#'
#'   n <- 100
#'   data1 <- GeneNet::ggm.simulate.data(n, pcor1)
#'   data2 <- GeneNet::ggm.simulate.data(n, pcor2)
#'
#'   X1 <- rbind(data1, data2) ## change-point at 100
#'   N1 <- detect.th(X1, approach = 'euclidean', scales = -1)
#'   N2 <- detect.th(X1, approach = 'infinity', scales = -1)
#'   N1$changepoints
#'   N1$time_series
#'   N2$changepoints
#'   N2$time_series
detect.th <- function(X, approach = c("euclidean", "infinity"),
                                        th_max = 2.25, th_sum = 0.65,
                                        pointsgen = 10, scales = -1,
                                        preaverage_gen = FALSE, scal_gen = 3,
                                        min_dist = 1) {
    approach = as.character(approach)
    sgn <- sign(stats::cor(X))
    if (!is.null(scales)) {
        stopifnot(sum(scales < 0) == length(scales))
        scales <- sort(scales, decreasing = TRUE)
    }
    if (ncol(X) <= 10) {
      approach[1] <- "infinity"
    }
    lsc <- length(scales)
    input <- t(hdbinseg::gen.input(t(X), scales, TRUE, diag = FALSE, sgn))
    if (approach[1] == "infinity") {
        cpt <- cpt_ts_Linf(input, thr_const_m = th_max,
                           points_m = pointsgen,
                           preproc = preaverage_gen, scl = scal_gen,
                           scale = scales)
    } else {
        cpt <- cpt_ts_L2(input, thr_const_s = th_sum,
                         points_s = pointsgen,
                         preproc = preaverage_gen, scl = scal_gen,
                         scale = scales)
    }
    no.of.cpt <- length(cpt$changepoints)
    d_cpt <- diff(cpt$changepoints)
    if ((min_dist <= 1) | (no.of.cpt <= 1) | (all(d_cpt > min_dist))) {
        return(list(changepoints = cpt$changepoints,
                    no.of.cpts = no.of.cpt,
                    time_series = cpt$time_series))
    } else {
        return(prune_cpt(input, cpts = cpt$changepoints,
                         distance = min_dist))
    }
}

#' Multiple change-point detection in the cross-covariance structure of
#' multivariate high-dimensional time series using a model selection
#' criterion optimisation
#'
#' This function detects multiple change-points in the cross-covariance
#' structure of a multivariate time series using a model selection
#' criterion optimisation.
#'
#' @export
#' @param X A numerical matrix representing the multivariate time series,
#'   with the columns representing its components.
#' @param approach A character string, which defines the metric to be used in
#'   order to detect the change-points. If approach = ``euclidean'', which is
#'   also the default value, then the \eqn{L_2} metric will be followed for
#'   the detection. If approach = ``infinity'', then the \eqn{L_{\infty}}
#'   metric will be used for the detection.
#' @param th_max A positive real number with default value equal to 2.1. It is
#'   used to define the threshold for the change-point overestimation step if
#'   the \eqn{L_{\infty}} metric is chosen in \code{approach} .
#' @param th_sum A positive real number with default value equal to 0.5. It is
#'   used to define the threshold for the change-point overestimation step if
#'   the \eqn{L_2} metric is chosen in \code{approach}.
#' @param pointsgen A positive integer with default value equal to 10. It
#'   defines the distance between two consecutive end- or start-points of
#'   the right- or left-expanding intervals, respectively; see Details
#'   for more information.
#' @param scales Negative integers for wavelet scales, with a small negative
#'   integer representing a fine scale. The default value is equal to -1.
#' @param alpha_gen A positive real number with default value equal to 0.1.
#'   It is used to define how strict the user wants to be with the penalty
#'   used.
#' @param preaverage_gen A logical variable with default value equal to
#'   \code{FALSE}. If \code{FALSE}, then pre-averaging the data is not
#'   required. If \code{TRUE}, then we need to pre-average the data before
#'   proceeding with the detection of the change-points.
#' @param scal_gen A positive integer number with default value equal to 3.
#'   It is used to define the way we pre-average the given data sequence
#'   only if \code{preaverage_gen = TRUE}. See the Details in
#'   \code{\link{preaverage}} for more information on how we pre-average.
#' @param min_dist A positive integer number with default value equal to 1. It
#'   is used in order to provide the minimum distance acceptable between
#'   detected change-points if such restrictions apply.
#' @details The time series \eqn{X_t} is of dimensionality \eqn{p} and we are
#'   looking for changes in the cross-covariance structure between the
#'   different time series components
#'   \eqn{X_{t}^{(1)}, X_{t}^{(2)}, ..., X_{t}^{(p)}}. We first use a
#'   wavelet-based approach for the various given scales in \code{scales} in
#'   order to transform the given time series \eqn{X_t} to a multiplicative
#'   model \eqn{Y_{t}^{(k)} = \sigma^{(k)}_t (Z_t^{(k)})^2; t=1,2,\ldots,T; k = 1,2,\ldots,d,}
#'   where \eqn{Z_t^{(k)}} is a sequence of standard normal random variables,
#'   \eqn{E(Y_t^{(k)}) = \sigma_t^{(k)}}, and \eqn{d} is the new
#'   dimensionality, which depends on the value given in \code{scales}.
#'   The function has been extensively tested on fMRI data, hence, its parameters
#'   have been tuned for this data type. The function might not work well in other
#'   structures, such as time series that are negatively serially correlated.
#' @return
#'   A list with the following components:
#'   \tabular{ll}{
#'    \cr \code{changepoints} \tab The locations of the detected change-points.
#'    \cr \code{no.of.cpts} \tab The number of the detected change-points.
#'    \cr \code{sol_path} \tab A vector containing the solution path.
#'    \cr \code{ic_curve}   \tab A vector with values of the information criterion for
#'    different number of change-points.
#'  }
#'   If the minimum distance between the detected change-points is less than
#'   the value given in \code{min_dist}, then only the number and the locations of the
#'   ``pruned'' change-points are returned.
#' @author Andreas Anastasiou, \email{anastasiou.andreas@ucy.ac.cy}
#' @references ``Cross-covariance isolate detect: a new change-point method
#'   for estimating dynamic functional connectivity'', Anastasiou et al (2020),
#'   preprint <doi:10.1101/2020.12.20.423696>.
#' @seealso \code{\link{detect.th}}.
#' @examples
#'   set.seed(11)
#'   A <- matrix(rnorm(10*200), nrow = 200) ## No change-point
#'   M1 <- detect.ic(A, approach = 'euclidean', scales = -1)
#'   M2 <- detect.ic(A, approach = 'infinity', scales = -1)
#'   M1$changepoints
#'   M2$changepoints
#'
#'   set.seed(1)
#'   num.nodes <- 30 # number of nodes
#'   etaA.1    <- 0.95
#'   etaA.2    <- 0.05
#'   pcor1     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.1)
#'   pcor2     <- GeneNet::ggm.simulate.pcor(num.nodes, etaA = etaA.2)
#'
#'   n <- 50
#'   data1 <- GeneNet::ggm.simulate.data(n, pcor1)
#'   data2 <- GeneNet::ggm.simulate.data(n, pcor2)

#'   X1 <- rbind(data1, data2, data1, data2) ## change-points at 50, 100, 150
#'   N1 <- detect.ic(X1, approach = 'euclidean', scales = -1)
#'   N2 <- detect.ic(X1, approach = 'infinity', scales = -1)
#'   N1$changepoints
#'   N2$changepoints
#'   N1$no.of.cpts
#'   N2$no.of.cpts
#'   N1$sol_path
#'   N2$sol_path
detect.ic <- function(X, approach = c("euclidean", "infinity"),
                                       th_max = 2.1, th_sum = 0.5,
                                       pointsgen = 10, scales = -1,
                                       alpha_gen = 0.1, preaverage_gen = FALSE,
                                       scal_gen = 3, min_dist = 1) {
  if (!(is.numeric(X))) {
    stop("The input in `X' should be a numeric matrix.")
  }
    sgn <- sign(stats::cor(X))
    if (!is.null(scales)) {
        stopifnot(sum(scales < 0) == length(scales))
        scales <- sort(scales, decreasing = TRUE)
    }
    if (ncol(X) <= 10) {
      approach[1] <- "infinity"
    }
    input <- t(hdbinseg::gen.input(t(X), scales, TRUE, diag = FALSE, sgn))
    input <- input + 10 ^ (-100)
    if (approach[1] == "infinity") {
        cpt <- ic_Linf(input, th_const = th_max,
                                   points = pointsgen, alpha = alpha_gen,
                                   preproc = preaverage_gen, scl = scal_gen)
    } else {
        cpt <- ic_L2(input, th_const = th_sum,
                                 points = pointsgen, alpha = alpha_gen,
                                 preproc = preaverage_gen, scl = scal_gen)
    }
    no.of.cpt <- cpt$no.of.cpts
    d_cpt <- diff(cpt$changepoints)
    if ((min_dist <= 1) | (no.of.cpt <= 1) | (all(d_cpt > min_dist))) {
        return(cpt)
    } else {
        return(prune_cpt(input, cpt$changepoints, distance = min_dist))
    }
}

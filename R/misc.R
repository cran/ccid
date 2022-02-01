options(expressions = 5e+05)
cusum <- function(x){
  if (!(is.numeric(x))) {
    stop("The input in `x' should be a numeric vector containing the data
         for which the CUSUM function will be calculated.")
  }
  n <- length(x)
  y <- cumsum(x)
  res <- sqrt(((n - 1):1)/n/(1:(n - 1))) * y[1:(n - 1)] - sqrt((1:(n - 1))/n/((n - 1):1)) * (y[n] - y[1:(n - 1)])
  return(res)
}

thr_Linf <- function(X, thr_const = 2.25,
                                 thr_fin = thr_const * sqrt(log(nrow(X))),
                                 s = 1, e = nrow(X), points = 10, k_l = 1,
                                 k_r = 1) {
  if (!(is.matrix(X))) {
    stop("The input in `X' should be a numeric matrix, with each data sequence
        we want to investigate being a column of this matrix.")
  }
  if ((thr_const <= 0) || (points <= 0)) {
    stop("The threshold constant as well as the `points' argument that
        represents the magnitude of the expansion for the intervals should
        be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps^0.5) {
    warning("The input for `points' should be a positive integer.
        If it is a positive real number then the integer part of the
        given number is used as the value of `points'.")
  }
  points <- as.integer(points)
  l <- length(X[, 1])
  r_e_points <- seq(points, l, points)
  l_e_points <- seq(l - points + 1, 1, -points)
  chp <- 0
  if (e - s <= 1) {
    cpt <- 0
  } else {
    pos_r <- numeric()
    CUSUM_r <- numeric()
    pos_l <- numeric()
    CUSUM_l <- numeric()
    moving_points <- IDetect::s_e_points(r_e_points, l_e_points, s, e)
    right_points <- moving_points[[1]]
    left_points <- moving_points[[2]]
    lur <- length(left_points)
    rur <- length(right_points)
    if (k_r < k_l) {
      while ((chp == 0) & (k_r < min(k_l, rur))) {
        lxr <- right_points[k_r] - s + 1
        ipcr <- matrix(NA, lxr - 1, ncol(X))
        mpcr <- numeric()
        for (i in seq_len(ncol(X))) {
          ipcr[, i] <- cusum(X[s:right_points[k_r], i]) / mean(X[s:right_points[k_r], i])
          mpcr[i] <- max(abs(ipcr[, i]))
        }
        pos.max <- which.max(mpcr)
        pos_r[k_r] <- which.max(abs(ipcr[, pos.max])) + s - 1
        CUSUM_r[k_r] <- mpcr[pos.max]
        if (CUSUM_r[k_r] > thr_fin) {
          chp <- pos_r[k_r]
        } else {
          k_r <- k_r + 1
        }
      }
    }
    if (k_l < k_r) {
      while ((chp == 0) & (k_l < min(k_r, lur))) {
        lxl <- e - left_points[k_l] + 1
        ipcl <- matrix(NA, lxl - 1, ncol(X))
        mpcl <- numeric()
        for (i in seq_len(ncol(X))) {
          ipcl[, i] <- cusum(X[left_points[k_l]:e, i]) / mean(X[left_points[k_l]:e, i])
          mpcl[i] <- max(abs(ipcl[, i]))
        }
        pos.max <- which.max(mpcl)
        pos_l[k_l] <- which.max(abs(ipcl[, pos.max])) + left_points[k_l] - 1
        CUSUM_l[k_l] <- mpcl[pos.max]
        if (CUSUM_l[k_l] > thr_fin) {
          chp <- pos_l[k_l]
        } else {
          k_l <- k_l + 1
        }
      }
    }
    if (chp == 0) {
      while ((chp == 0) & (k_l <= lur) & (k_r <= rur)) {
        lxr <- right_points[k_r] - s + 1
        ipcr <- matrix(NA, lxr - 1, ncol(X))
        mpcr <- numeric()
        for (i in seq_len(ncol(X))) {
          ipcr[, i] <- cusum(X[s:right_points[k_r], i]) / mean(X[s:right_points[k_r], i])
          mpcr[i] <- max(abs(ipcr[, i]))
        }
        pos.max <- which.max(mpcr)
        pos_r[k_r] <- which.max(abs(ipcr[, pos.max])) + s - 1
        CUSUM_r[k_r] <- mpcr[pos.max]
        if (CUSUM_r[k_r] > thr_fin) {
          chp <- pos_r[k_r]
        } else {
          lxl <- e - left_points[k_l] + 1
          ipcl <- matrix(NA, lxl - 1, ncol(X))
          mpcl <- numeric()
          for (i in seq_len(ncol(X))) {
            ipcl[, i] <- cusum(X[left_points[k_l]:e, i]) / mean(X[left_points[k_l]:e, i])
            mpcl[i] <- max(abs(ipcl[, i]))
          }
          pos.max <- which.max(mpcl)
          pos_l[k_l] <- which.max(abs(ipcl[, pos.max])) + left_points[k_l] - 1
          CUSUM_l[k_l] <- mpcl[pos.max]
          if (CUSUM_l[k_l] > thr_fin) {
            chp <- pos_l[k_l]
          } else {
            k_r <- k_r + 1
            k_l <- k_l + 1
          }
        }
      }
    }
    if (chp != 0) {
      if (chp > ((e + s) / 2)) {
        r <- thr_Linf(X, s = s, e = chp, points = points,
                                  thr_fin = thr_fin, k_r = k_r,
                                  k_l = 1)
      } else {
        r <- thr_Linf(X, s = chp + 1, e = e,
                                  points = points, thr_fin = thr_fin,
                                  k_r = 1, k_l = max(1, k_l - 1))
      }
      cpt <- c(chp, r)
    } else {
      cpt <- chp
    }
  }
  cpt <- cpt[cpt != 0]
  return(sort(cpt))
}


thr_L2 <- function(X, thr_const = 0.65,
                               thr_fin = thr_const * sqrt(log(nrow(X))),
                               s = 1, e = nrow(X), points = 10, k_l = 1,
                               k_r = 1) {
  if (!(is.matrix(X))) {
    stop("The input in `X' should be a numeric matrix, with each data
        sequence we want to investigate being a column of this matrix.")
  }
  if ((thr_fin <= 0) || (points <= 0)) {
    stop("The threshold constant as well as the `points' argument that
        represents the magnitude of the expansion for the intervals should be
        positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps^0.5) {
    warning("The input for `points' should be a positive integer. If it is
        a positive real number then the integer part of the given number is
        used as the value of `points'.")
  }
  points <- as.integer(points)
  l <- length(X[, 1])
  r_e_points <- seq(points, l, points)
  l_e_points <- seq(l - points + 1, 1, -points)
  chp <- 0
  if (e - s <= 1) {
    cpt <- 0
  } else {
    pos_r <- numeric()
    CUSUM_r <- numeric()
    pos_l <- numeric()
    CUSUM_l <- numeric()
    moving_points <- IDetect::s_e_points(r_e_points, l_e_points, s, e)
    right_points <- moving_points[[1]]
    left_points <- moving_points[[2]]
    lur <- length(left_points)
    rur <- length(right_points)
    if (k_r < k_l) {
      while ((chp == 0) & (k_r < min(k_l, rur))) {
        lxr <- right_points[k_r] - s + 1
        ipcr <- matrix(NA, lxr - 1, ncol(X))
        for (i in seq_len(ncol(X))) {
          ipcr[, i] <- cusum(X[s:right_points[k_r], i]) / mean(X[s:right_points[k_r], i])
        }
        sipcr <- sqrt(rowSums(ipcr^2)) / sqrt(ncol(X))
        pos_r[k_r] <- which.max(sipcr) + s - 1
        CUSUM_r[k_r] <- sipcr[pos_r[k_r] - s + 1]
        if (CUSUM_r[k_r] > thr_fin) {
          chp <- pos_r[k_r]
        } else {
          k_r <- k_r + 1
        }
      }
    }
    if (k_l < k_r) {
      while ((chp == 0) & (k_l < min(k_r, lur))) {
        lxl <- e - left_points[k_l] + 1
        ipcl <- matrix(NA, lxl - 1, ncol(X))
        for (i in seq_len(ncol(X))) {
          ipcl[, i] <- cusum(X[left_points[k_l]:e, i]) / mean(X[left_points[k_l]:e, i])
        }
        sipcl <- sqrt(rowSums(ipcl^2)) / sqrt(ncol(X))
        pos_l[k_l] <- which.max(sipcl) + left_points[k_l] - 1
        CUSUM_l[k_l] <- sipcl[pos_l[k_l] - left_points[k_l] + 1]
        if (CUSUM_l[k_l] > thr_fin) {
          chp <- pos_l[k_l]
        } else {
          k_l <- k_l + 1
        }
      }
    }
    if (chp == 0) {
      while ((chp == 0) & (k_l <= lur) & (k_r <= rur)) {
        lxr <- right_points[k_r] - s + 1
        ipcr <- matrix(NA, lxr - 1, ncol(X))
        for (i in seq_len(ncol(X))) {
          ipcr[, i] <- cusum(X[s:right_points[k_r], i]) / mean(X[s:right_points[k_r], i])
        }
        sipcr <- sqrt(rowSums(ipcr^2)) / sqrt(ncol(X))
        pos_r[k_r] <- which.max(sipcr) + s - 1
        CUSUM_r[k_r] <- sipcr[pos_r[k_r] - s + 1]
        if (CUSUM_r[k_r] > thr_fin) {
          chp <- pos_r[k_r]
        } else {
          lxl <- e - left_points[k_l] + 1
          ipcl <- matrix(NA, lxl - 1, ncol(X))
          for (i in seq_len(ncol(X))) {
            ipcl[, i] <- cusum(X[left_points[k_l]:e, i]) / mean(X[left_points[k_l]:e, i])
          }
          sipcl <- sqrt(rowSums(ipcl^2)) / sqrt(ncol(X))
          pos_l[k_l] <- which.max(sipcl) + left_points[k_l] - 1
          CUSUM_l[k_l] <- sipcl[pos_l[k_l] - left_points[k_l] + 1]
          if (CUSUM_l[k_l] > thr_fin) {
            chp <- pos_l[k_l]
          } else {
            k_r <- k_r + 1
            k_l <- k_l + 1
          }
        }
      }
    }
    if (chp != 0) {
      if (chp > ((e + s) / 2)) {
        r <- thr_L2(X, s = s, e = chp, points = points,
                                thr_fin = thr_fin, k_r = k_r,
                                k_l = 1)
      } else {
        r <- thr_L2(X, s = chp + 1, e = e, points = points,
                                thr_fin = thr_fin, k_r = 1,
                                k_l = max(1, k_l - 1))
      }
      cpt <- c(chp, r)
    } else {
      cpt <- chp
    }
  }
  cpt <- cpt[cpt != 0]
  return(sort(cpt))
}

cpt_ts_Linf <- function(X_m, thr_const_m = 2.25, points_m = 10,
                              preproc = FALSE, scl = 3, scale = -1) {
  if (preproc) {
    X_m <- preaverage(X_m, scal = scl)
  }
  cpt_m <- thr_Linf(X = X_m, thr_const = thr_const_m,
                                points = points_m)
  ts_m <- match.cpt.ts(X = X_m, cpt = cpt_m, scales = scale)
  if (preproc) {
    cpt_m <- (cpt_m - 1) * scl + round(scl / 2)
}
  return(list(changepoints = cpt_m, time_series = ts_m))
}

cpt_ts_L2 <- function(X_s, thr_const_s = 0.65, points_s = 10,
                            preproc = FALSE, scl = 3, scale = -1) {
  if (preproc) {
    X_s <- preaverage(X_s, scal = scl)
  }
  cpt_s <- thr_L2(X = X_s, thr_const = thr_const_s,
                              points = points_s)
  ts_s <- match.cpt.ts(X = X_s, cpt = cpt_s, scales = scale)
  if (preproc) {
    cpt_s <- (cpt_s - 1) * scl + round(scl / 2)
  }
  return(list(changepoints = cpt_s, time_series = ts_s))
}

sol_path_Linf <- function(X, thr_ic = 2.1, points = 10) {
  if (!(is.numeric(X))) {
    stop("The input in `X' should be a numeric matrix containing the
        multivariate data for which the solution path will be calculated.")
  }
  if ((thr_ic <= 0) || (points <= 0)) {
    stop("The threshold constant as well as the `points' argument that
        represents the magnitude of the expansion for the intervals should
        be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps^0.5) {
    warning("The input for `points' should be a positive integer. If it
        is a positive real number then the integer part of the given number
        is used as the value of `points'.")
  }
  lx_ic <- nrow(X)
  points <- as.integer(points)
  cpt_lower <- thr_Linf(X, thr_const = thr_ic, points = points)
  lcpt_ic <- length(cpt_lower)
  cs <- cumul_sum_matrix(X)
  if ((lcpt_ic == 1) | (lcpt_ic == 0)) {
    return(list(solution_path = cpt_lower))
  } else {
    seb_set <- c(unique(c(1, cpt_lower)), lx_ic)
    lseb_set <- length(seb_set)
    min_C <- numeric()
    while (lseb_set >= 3) {
      Rs <- cusum_point(X, y = cs, seb_set[1:(lseb_set - 2)],
                              seb_set[3:(lseb_set)],
                              seb_set[2:(lseb_set - 1)])
      min_Rs <- seb_set[2:(lseb_set - 1)][which.min(Rs)]
      min_C <- c(min_C, min_Rs)
      seb_set <- seb_set[-which(seb_set == min_Rs)]
      lseb_set <- lseb_set - 1
    }
    return(list(solution_path = min_C[rev(seq_len(length(min_C)))]))
  }
}


sol_path_L2 <- function(X, thr_ic = 0.5, points = 10) {
  if (!(is.numeric(X))) {
    stop("The input in `X' should be a numeric matrix containing the
        multivariate data for which the solution path will be calculated.")
  }
  if ((thr_ic <= 0) || (points <= 0)) {
    stop("The threshold constant as well as the `points' argument that
        represents the magnitude of the expansion for the intervals should
        be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps^0.5) {
    warning("The input for `points' should be a positive integer. If it is
        a positive real number then the integer part of the given number is
        used as the value of `points'.")
  }
  lx_ic <- nrow(X)
  points <- as.integer(points)
  cpt_lower <- thr_L2(X, thr_const = thr_ic, points = points)
  lcpt_ic <- length(cpt_lower)
  cs <- cumul_sum_matrix(X)
  if ((lcpt_ic == 1) | (lcpt_ic == 0)) {
    return(list(solution_path = cpt_lower))
  } else {
    seb_set <- c(unique(c(1, cpt_lower)), lx_ic)
    lseb_set <- length(seb_set)
    min_C <- numeric()
    while (lseb_set >= 3) {
      Rs <- cusum_point(X, y = cs, seb_set[1:(lseb_set - 2)],
                              seb_set[3:(lseb_set)],
                              seb_set[2:(lseb_set - 1)])
      min_Rs <- seb_set[2:(lseb_set - 1)][which.min(Rs)]
      min_C <- c(min_C, min_Rs)
      seb_set <- seb_set[-which(seb_set == min_Rs)]
      lseb_set <- lseb_set - 1
    }
    sp <- min_C[length(min_C):1]
    return(list(solution_path = sp))
  }
}

ic_Linf <- function(X, th_const = 2.1, Kmax = 100, points = 10,
                                alpha = 0.1, preproc = FALSE, scl = 3) {
  if (!(is.numeric(X))) {
    stop("The input in `X' should be a numeric matrix containing the
        multivariate data in which you want to look for change-points.")
  }
  if ((th_const <= 0) || (points <= 0)) {
    stop("The threshold constant as well as the `points' argument that
        represents the magnitude of the expansion for the intervals should
        be positive numbers.")
  }
  if (preproc) {
    X <- preaverage(X, scal = scl)
  }
  result <- list()
  lx <- nrow(X)
  ld <- ncol(X)
  if (Kmax == 0 || lx <= 3) {
    result$changepoints <- NA
    result$no.of.cpts <- 0
    result$sol_path <- NA
    result$ic_curve <- NA
    if (Kmax == 0) {
      stop("No change-points found, choose larger Kmax")
    } else {
      stop("Sample size is too small")
    }
  } else {
    cpt_cand <- sol_path_Linf(X, thr_ic = th_const, points = points)[[1]]
    if (length(cpt_cand) == 0) {
      result$changepoints <- NA
      result$no.of.cpts <- as.integer(0)
      result$sol_path <- NA
      result$ic_curve <- NA
    }
    else{
      if (length(cpt_cand) > min(Kmax, lx - 2)) {
        cpt_cand <- cpt_cand[1:min(Kmax, lx - 2)]
      }
      len_cpt <- length(cpt_cand)
      result$changepoints <- numeric()
      result$no.of.cpts <- integer()
      result$sol_path <- cpt_cand
      result$ic_curve <- rep(0, len_cpt + 1)
      if (len_cpt) {
        for (i in len_cpt:1) {
          min_log_lik <- loglik_chisq(X, cpt = cpt_cand[1:i])
          result$ic_curve[i + 1] <- min_log_lik + ((i + 1) * ld * (log(lx)) ^ (alpha)) / 2
        }
      }
      result$ic_curve[1] <- loglik_chisq(X, cpt = NULL) +  (ld * (log(lx)) ^ (alpha)) / 2
      tmp <- stats::quantile(which.min(result$ic_curve), 0.5, type = 3)
      if (tmp == 1) {
        result$changepoints <- NA
        result$no.of.cpts <- as.integer(0)
      } else {
        result$changepoints <- sort(cpt_cand[1:(tmp - 1)])
        result$no.of.cpts <- as.integer(tmp - 1)
      }
      if (preproc) {
        result$changepoints <- (result$changepoints - 1) * scl + round(scl / 2)
        result$sol_path <- (result$sol_path - 1) * scl + round(scl / 2)
      }
    }
    return(result)
  }
}

ic_L2 <- function(X, th_const = 0.5, Kmax = 100, points = 10,
                              alpha = 0.1, preproc = FALSE, scl = 3) {

  if (!(is.numeric(X))) {
    stop("The input in `X' should be a numeric matrix containing the
        multivariate data in which you want to look for change-points.")
  }
  if ((th_const <= 0) || (points <= 0)) {
    stop("The threshold constant as well as the `points' argument that
        represents the magnitude of the expansion for the intervals should
        be positive numbers.")
  }
  if (preproc) {
    X <- preaverage(X, scal = scl)
  }
  result <- list()
  lx <- nrow(X)
  ld <- ncol(X)
  if (Kmax == 0 || lx <= 3) {
    result$changepoints <- NA
    result$no.of.cpts <- 0
    result$sol_path <- NA
    result$ic_curve <- NA
    if (Kmax == 0) {
      stop("No change-points found, choose larger Kmax")
    } else {
      stop("Sample size is too small")
    }
  } else {
    cpt_cand <- sol_path_L2(X, thr_ic = th_const, points = points)[[1]]
    if (length(cpt_cand) == 0) {
      result$changepoints <- NA
      result$no.of.cpts <- as.integer(0)
      result$sol_path <- NA
      result$ic_curve <- NA
    } else {
      if (length(cpt_cand) > min(Kmax, lx - 2)) {
        cpt_cand <- cpt_cand[1:min(Kmax, lx - 2)]
      }
      result$changepoints <- numeric()
      result$no.of.cpts <- integer()
      result$sol_path <- cpt_cand
      len_cpt <- length(cpt_cand)
      result$ic_curve <- rep(0, len_cpt + 1)
      for (i in len_cpt:1) {
        min_log_lik <- loglik_chisq(X, cpt = cpt_cand[1:i])
        result$ic_curve[i + 1] <- min_log_lik + ((i + 1) * ld * (log(lx)) ^ (alpha)) / 2
      }
      result$ic_curve[1] <- loglik_chisq(X, cpt = NULL) + (ld * (log(lx)) ^ (alpha)) / 2
      tmp <- stats::quantile(which.min(result$ic_curve), 0.5, type = 3)
      if (tmp == 1) {
        result$changepoints <- NA
        result$no.of.cpts <- as.integer(0)
      } else {
        result$changepoints <- sort(cpt_cand[1:(tmp - 1)])
        result$no.of.cpts <- as.integer(tmp - 1)
      }
      if (preproc) {
        result$changepoints <- (result$changepoints - 1) * scl + round(scl / 2)
        result$sol_path <- (result$sol_path - 1) * scl + round(scl / 2)
      }
    }
    return(result)
  }
}

cusum_one_multi <- function(x, s, e, b) {
    if (!(is.numeric(x))) {
        stop("The input in `x' should be a numeric vector.")
    }
    y <- cumsum(x)
    l <- numeric()
    d <- numeric()
    result <- numeric()
    if ((length(s) != length(b)) || (length(s) != length(e)) || (length(e) != length(b))) {
        stop("The vectors s, b, e, should be of the same length")
    }
    if (any(s < 1) | any(b < 1) | any(e < 1)) {
        stop("The entries of the vectors s, b, e should be positive integers.")
    }
    if (any(s > b) | any(b >= e)) {
        stop("The value for b should be in the interval [s,e)")
    }
    if ((any(abs((s - round(s))) > .Machine$double.eps^0.5)) || (any(abs((b - round(b))) > .Machine$double.eps^0.5)) || (any(abs((e -
        round(e))) > .Machine$double.eps^0.5))) {
        stop("The input values  for s, b, and  e should be positive integers.")
    }
    for (j in seq_len(length(b))) {
        l[j] <- e[j] - s[j] + 1
        d[j] <- ifelse(s[j] == 1, 0, y[s[j] - 1])
        result[j] <- (1 / mean(x[s[j]:e[j]])) * (abs(sqrt((e[j] - b[j]) / (l[j] * (b[j] - s[j] + 1))) * (y[b[j]] - d[j]) - sqrt((b[j] - s[j] + 1) / (l[j] * (e[j] - b[j]))) * (y[e[j]] - y[b[j]])))
    }
    return(result)
}

cumul_sum_matrix <- function(x) {
    lx <- dim(x)[1]
    l_grid <- dim(x)[2]
    cumul_sum <- matrix(NA, lx, l_grid)
    i <- 1
    while (i <= l_grid) {
        cumul_sum[, i] <- cumsum(x[, i])
        i <- i + 1
    }
    return(cumul_sum)
}

cusum_point <- function(X, y = cumul_sum_matrix(X), s, e, b) {
    if (!(is.numeric(X))) {
        stop("The input in `X' should be a numeric matrix")
    }
    if ((length(s) != length(b)) || (length(s) != length(e)) || (length(e) != length(b))) {
        stop("The vectors s, b, e, should be of the same length")
    }
    if (any(s < 1) | any(b < 1) | any(e < 1)) {
        stop("The entries of the vectors s, b, e should be positive integers.")
    }
    if (any(s > b) | any(b >= e)) {
        stop("The value for b should be in the interval [s,e)")
    }
    if ((any(abs((s - round(s))) > .Machine$double.eps^0.5)) || (any(abs((b - round(b))) > .Machine$double.eps^0.5)) || (any(abs((e -
        round(e))) > .Machine$double.eps^0.5))) {
        stop("The input values  for s, b, and  e should be positive integers.")
    }
    l_grid <- dim(y)[2]
    l <- numeric()
    alpha <- numeric()
    beta <- numeric()
    result_matrix <- matrix(NA, length(b), l_grid)
    d <- matrix(NA, length(b), l_grid)
    for (j in seq_len(length(b))) {
        l[j] <- e[j] - s[j] + 1
        alpha[j] <- b[j] - s[j] + 1
        beta[j] <- e[j] - b[j]
        if (s[j] == 1) {
            d[j, ] <- rep(0, l_grid)
        } else {
            d[j, ] <- y[s[j] - 1, ]
        }
        for (i in 1:l_grid) {
            result_matrix[j, i] <- (1 / mean(X[s[j]:e[j], i])) * abs(sqrt(beta[j] / (l[j] * alpha[j])) * (y[b[j], i] - d[j, i]) - sqrt(alpha[j] / (l[j] * beta[j])) * (y[e[j], i] - y[b[j], i]))
        }
    }
    return(apply(result_matrix, 1, max))
}

loglik_chisq <- function(X, cpt, no.cpt = length(cpt)) {
    lx <- nrow(X)
    ld <- ncol(X)
    result_temp <- numeric()
    temp_log <- numeric()
    if (no.cpt == 0) {
        for (j in 1:ld) {
            result_temp[j] <- lx / 2 + (lx / 2) * log(pi * mean(X[, j])) + 0.5 * sum(log(2 * X[, j]))
        }
        result <- sum(result_temp)
    } else {
        crucial_point <- c(0, sort(cpt), lx)
        dcpt <- diff(crucial_point)
        log_s <- numeric()
        for (j in 1:ld) {
            for (k in 1:(no.cpt + 1)) {
                temp_log[k] <- (dcpt[k] / 2) * log(pi * mean(X[(crucial_point[k] + 1):crucial_point[k + 1], j])) + 0.5 * sum(log(2 * X[(crucial_point[k] + 1):crucial_point[k + 1], j])) + dcpt[k] / 2
            }
            log_s[j] <- sum(temp_log)
        }
        result <- sum(log_s)
    }
    result
}


prune_cpt <- function(X, cpts, distance) {
    if (cpts[1] == 1){
      cpts <- cpts[-1]
    }
  if (cpts[length(cpts)] == nrow(X)){
    cpts <- cpts[-length(cpts)]
  }
    l_cpts <- length(cpts)
    dc <- diff(cpts)
    ind <- which(dc < distance)
    l_ind <- length(ind)
    if (l_ind == 0) {
        return(list(changepoints = cpts, no.of.cpts = l_cpts))
    } else {
        cands <- sort(cpts[unique(c(ind, ind + 1))])
        l_cands <- length(cands)
        cs <- cumul_sum_matrix(X)
        seb_set <- numeric()
        seb_set[1] <- ifelse(cands[1] == cpts[1], 1, cpts[ind[1] - 1])
        seb_set[l_cands + 2] <- ifelse(cands[l_cands] == cpts[l_cpts], nrow(X), cpts[ind[l_ind] + 2])
        seb_set <- c(seb_set[1], cands, seb_set[l_cands + 2])
        lseb_set <- length(seb_set)
        min_C <- numeric()
        min_Rs <- numeric()
        while (l_ind >= 1) {
            Rs <- cusum_point(X, y = cs, seb_set[1:(lseb_set - 2)], seb_set[3:(lseb_set)], seb_set[2:(lseb_set - 1)])
            min_Rs <- cands[which.min(Rs)]
            min_C <- c(min_C, min_Rs)
            seb_set <- seb_set[-which(seb_set == min_Rs)]
            lseb_set <- length(seb_set)
            cands <- seb_set[2:(lseb_set - 1)]
            dc <- diff(cands)
            ind <- which(dc < distance)
            l_ind <- length(ind)
            cands <- sort(cands[unique(c(ind, ind + 1))])
            l_cands <- length(cands)
            indicator <- which(seb_set %in% cands)
            seb_set <- c(seb_set[indicator[1] - 1], cands, seb_set[indicator[length(indicator)] + 1])
            lseb_set <- length(seb_set)
          }
        cpts <- cpts[!cpts %in% min_C]
        l.cpts <- length(cpts)
        return(list(changepoints = cpts, no.of.cpts = l.cpts))
    }
}


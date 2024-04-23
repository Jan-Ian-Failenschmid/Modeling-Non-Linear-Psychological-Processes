#' ----------------------------------------------------------------------------#
#' Title: Local polynomial methods                                             #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 11-04-2024                                                    #
#' -----                                                                       #
#' Last Modified: 23-04-2024                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#

#' Class and methods for the local polynomial method
require(nprobust)

# Locpol regression
setClass(
  "method_locpol",
  contains = "method"
)

setMethod("fit", "method_locpol", function(method, data) {
  loc_fit_list <- lapply(c(1, 3, 5), function(p, method, data) {
    loc_fit <- lprobust_cust(
      x = data$time, y = data$y_obs, eval = data$time,
      p = p, kernel = "gau", bwselect = "imse-dpi",
      bwcheck = 0, diag_A = TRUE
    )
    gcv <- get_cv(loc_fit$Estimate, data$y_obs)
    return(list(loc_fit, gcv))
  }, method = method, data = data)

  loc_fit <- loc_fit_list[[which.min(sapply(loc_fit_list, "[[", 2))]][[1]]

  # Method generics schould always return the aFdjusted method object
  slot(method, "fit") <- loc_fit
  slot(method, "converged") <- TRUE

  return(method)
})

setMethod(
  "calculate_performance_measures", "method_locpol",
  function(method, data) {
    if (slot(method, "converged")) {
      inf <- as.data.frame(slot(method, "fit")$Estimate)

      # Get mean inference
      slot(method, "estimate") <- inf$tau.bc
      slot(method, "ci") <- list(
        ub = as.vector(inf$tau.bc + (qnorm(0.975) * inf$se.rb)),
        lb = as.vector(inf$tau.bc - (qnorm(0.975) * inf$se.rb))
      )

      # Calculate mse
      slot(method, "mse") <- calc_mse(method, data)

      # Calculate gcv
      slot(method, "gcv") <- get_cv(
        slot(method, "fit")$Estimate,
        data$y_obs
      )

      # Calculate confidence interval coverage
      slot(method, "ci_coverage") <- ci_test(method, data)
    } else {
      method <- set_na(method)
    }

    # Method generics schould always return the adjusted method object
    return(method)
  }
)

lprobust_cust <- function(
    y, x, eval = NULL, neval = NULL, p = NULL, deriv = NULL,
    h = NULL, b = NULL, rho = 1, kernel = "epa", bwselect = NULL,
    bwcheck = 21, bwregul = 1, imsegrid = 30, vce = "nn", covgrid = FALSE,
    cluster = NULL, nnmatch = 3, level = 95, interior = FALSE,
    subset = NULL, diag_A = TRUE) {
  #' Customization of the lprobust function from the nprobust package. 
  #' Code added to obtain analytic cross-validatin criteria during the 
  #' estimation. 
  if (!is.null(subset)) {
    x <- x[subset]
    y <- y[subset]
  }
  na.ok <- complete.cases(x) & complete.cases(y)
  if (!is.null(cluster)) {
    if (!is.null(subset)) {
      cluster <- cluster[subset]
    }
    na.ok <- na.ok & complete.cases(cluster)
  }
  x <- x[na.ok]
  y <- y[na.ok]
  if (!is.null(cluster)) {
    cluster <- cluster[na.ok]
  }
  if (!is.null(deriv) & is.null(p)) {
    p <- deriv + 1
  }
  if (is.null(p)) {
    p <- 1
  }
  if (is.null(deriv)) {
    deriv <- 0
  }
  q <- p + 1
  x.max <- max(x)
  x.min <- min(x)
  N <- length(x)
  if (is.null(eval)) {
    if (is.null(neval)) {
      eval <- seq(x.min, x.max, length.out = 30)
    } else {
      eval <- seq(x.min, x.max, length.out = neval)
    }
  }
  neval <- length(eval)
  if (is.null(h) & is.null(bwselect) & neval == 1) {
    bwselect <- "mse-dpi"
  }
  if (is.null(h) & is.null(bwselect) & neval > 1) {
    bwselect <- "imse-dpi"
  }
  if (vce == "nn") {
    order.x <- order(x)
    x <- x[order.x]
    y <- y[order.x]
    if (!is.null(cluster)) {
      cluster <- cluster[order.x]
    }
  }
  kernel <- tolower(kernel)
  bwselect <- tolower(bwselect)
  vce <- tolower(vce)
  vce_type <- "NN"
  if (vce == "hc0") {
    vce_type <- "HC0"
  }
  if (vce == "hc1") {
    vce_type <- "HC1"
  }
  if (vce == "hc2") {
    vce_type <- "HC2"
  }
  if (vce == "hc3") {
    vce_type <- "HC3"
  }
  if (vce == "cluster") {
    vce_type <- "Cluster"
  }
  if (vce == "nncluster") {
    vce_type <- "NNcluster"
  }
  exit <- 0
  if (kernel != "gau" & kernel != "gaussian" & kernel != "uni" &
    kernel != "uniform" & kernel != "tri" & kernel != "triangular" &
    kernel != "epa" & kernel != "epanechnikov" & kernel !=
    "") {
    print("kernel incorrectly specified")
    exit <- 1
  }
  if (vce != "nn" & vce != "" & vce != "hc1" & vce != "hc2" &
    vce != "hc3" & vce != "hc0") {
    print("vce incorrectly specified")
    exit <- 1
  }
  if (p < 0 | deriv < 0 | nnmatch <= 0) {
    print("p,q,deriv and matches should be positive integers")
    exit <- 1
  }
  if (deriv > p) {
    print("deriv can only be equal or lower p")
    exit <- 1
  }
  if (level > 100 | level <= 0) {
    print("level should be set between 0 and 100")
    exit <- 1
  }
  if (rho < 0) {
    print("rho should be greater than 0")
    exit <- 1
  }
  if (exit > 0) {
    stop()
  }
  if (!is.null(h)) {
    bwselect <- "Manual"
  }
  if (!is.null(h) & rho > 0 & is.null(b)) {
    b <- h / rho
  }
  kernel.type <- "Epanechnikov"
  if (kernel == "triangular" | kernel == "tri") {
    kernel.type <- "Triangular"
  }
  if (kernel == "uniform" | kernel == "uni") {
    kernel.type <- "Uniform"
  }
  if (kernel == "gaussian" | kernel == "gau") {
    kernel.type <- "Gaussian"
  }
  if (is.null(h)) {
    lpbws <- lpbwselect(
      y = y, x = x, eval = eval, deriv = deriv,
      p = p, vce = vce, cluster = cluster, bwselect = bwselect,
      interior = interior, kernel = kernel, bwcheck = bwcheck,
      bwregul = bwregul, imsegrid = imsegrid, subset = subset
    )
    h <- lpbws$bws[, 2]
    b <- lpbws$bws[, 3]
    if (rho > 0) {
      b <- h / rho
    }
    rho <- h / b
  }
  if (length(h) == 1 & neval > 1) {
    h <- rep(h, neval)
    b <- rep(b, neval)
    rho <- h / b
  }
  dups <- dupsid <- 0
  if (vce == "nn") {
    for (j in 1:N) {
      dups[j] <- sum(x == x[j])
    }
    j <- 1
    while (j <= N) {
      dupsid[j:(j + dups[j] - 1)] <- 1:dups[j]
      j <- j + dups[j]
    }
  }
  cov.p <- NULL
  if (!diag_A) {
    Estimate <- matrix(NA, neval, 8)
    colnames(Estimate) <- c(
      "eval", "h", "b", "N", "tau.us", "tau.bc",
      "se.us", "se.rb"
    )
  } else {
    Estimate <- matrix(NA, neval, 10)
    colnames(Estimate) <- c(
      "eval", "h", "b", "N", "tau.us", "tau.bc",
      "se.us", "se.rb", "A.us", "A.bc"
    )
  }

  for (i in 1:neval) {
    if (!is.null(bwcheck)) {
      bw.min <- sort(abs(x - eval[i]))[bwcheck]
      h[i] <- max(h[i], bw.min)
      b[i] <- max(b[i], bw.min)
    }
    w.h <- W.fun((x - eval[i]) / h[i], kernel) / h[i]
    w.b <- W.fun((x - eval[i]) / b[i], kernel) / b[i]
    ind.h <- w.h > 0
    ind.b <- w.b > 0
    N.h <- sum(ind.h)
    N.b <- sum(ind.b)
    ind <- ind.b
    if (h[i] > b[i]) {
      ind <- ind.h
    }
    if (diag_A) {
      A_ind <- numeric(length(y))
      A_ind[i] <- 1
      eA_ind <- A_ind[ind]
    }
    eN <- sum(ind)
    eY <- y[ind]
    eX <- x[ind]
    W.h <- w.h[ind]
    W.b <- w.b[ind]
    eC <- NULL
    if (!is.null(cluster)) {
      eC <- cluster[ind]
    }
    edups <- edupsid <- 0
    if (vce == "nn") {
      edups <- dups[ind]
      edupsid <- dupsid[ind]
    }
    u <- (eX - eval[i]) / h[i]
    R.q <- matrix(NA, eN, (q + 1))
    for (j in 1:(q + 1)) {
      R.q[, j] <- (eX - eval[i])^(j -
        1)
    }
    R.p <- R.q[, 1:(p + 1)]
    L <- crossprod(R.p * W.h, u^(p + 1))
    invG.q <- qrXXinv((sqrt(W.b) * R.q))
    invG.p <- qrXXinv((sqrt(W.h) * R.p))
    e.p1 <- matrix(0, (q + 1), 1)
    e.p1[p + 2] <- 1
    e.v <- matrix(0, (p + 1), 1)
    e.v[deriv + 1] <- 1
    Q.q <- t(t(R.p * W.h) - h[i]^(p + 1) * (L %*% t(e.p1)) %*%
      t(t(invG.q %*% t(R.q)) * W.b))
    beta.p <- invG.p %*% crossprod(R.p * W.h, eY)
    beta.q <- invG.q %*% crossprod(R.q * W.b, eY)
    beta.bc <- invG.p %*% crossprod(Q.q, eY)
    tau.cl <- factorial(deriv) * beta.p[(deriv + 1), 1]
    tau.bc <- factorial(deriv) * beta.bc[(deriv + 1), 1]
    if (diag_A) {
      A.cl <- factorial(deriv) *
        (invG.p %*% crossprod(R.p * W.h, eA_ind))[(deriv + 1), 1]
      A.bc <- factorial(deriv) *
        (invG.p %*% crossprod(Q.q, eA_ind))[(deriv + 1), 1]
    }
    hii <- predicts.p <- predicts.q <- 0
    if (vce == "hc0" | vce == "hc1" | vce == "hc2" | vce ==
      "hc3") {
      predicts.p <- R.p %*% beta.p
      predicts.q <- R.q %*% beta.q
      if (vce == "hc2" | vce == "hc3") {
        hii <- matrix(NA, eN, 1)
        for (j in 1:eN) {
          hii[j] <- R.p[j, ] %*% invG.p %*%
            (R.p * W.h)[j, ]
        }
      }
    }
    res.h <- lprobust.res(
      eX, eY, predicts.p, hii, vce, nnmatch,
      edups, edupsid, p + 1
    )
    if (vce == "nn") {
      res.b <- res.h
    } else {
      res.b <- lprobust.res(
        eX, eY, predicts.q, hii, vce,
        nnmatch, edups, edupsid, q + 1
      )
    }
    V.Y.cl <- invG.p %*% lprobust.vce(
      as.matrix(R.p * W.h),
      res.h, eC
    ) %*% invG.p
    V.Y.bc <- invG.p %*% lprobust.vce(Q.q, res.b, eC) %*%
      invG.p
    se.cl <- sqrt(factorial(deriv)^2 * V.Y.cl[
      deriv + 1,
      deriv + 1
    ])
    se.rb <- sqrt(factorial(deriv)^2 * V.Y.bc[
      deriv + 1,
      deriv + 1
    ])
    if (!diag_A) {
      Estimate[i, ] <- c(
        eval[i], h[i], b[i], eN, tau.cl, tau.bc,
        se.cl, se.rb
      )
    } else {
      Estimate[i, ] <- c(
        eval[i], h[i], b[i], eN, tau.cl, tau.bc,
        se.cl, se.rb, A.cl, A.bc
      )
    }
  }
  cov.us <- cov.rb <- matrix(NA, neval, neval)
  if (covgrid == TRUE) {
    for (i in 1:neval) {
      for (j in i:neval) {
        w.h.i <- W.fun((x - eval[i]) / h[i], kernel) / h[i]
        w.b.i <- W.fun((x - eval[i]) / b[i], kernel) / b[i]
        ind.h.i <- w.h.i > 0
        ind.b.i <- w.b.i > 0
        N.h.i <- sum(ind.h.i)
        N.b.i <- sum(ind.b.i)
        ind.i <- ind.b.i
        if (h[i] > b[i]) {
          ind.i <- ind.h.i
        }
        w.h.j <- W.fun((x - eval[j]) / h[j], kernel) / h[j]
        w.b.j <- W.fun((x - eval[j]) / b[j], kernel) / b[j]
        ind.h.j <- w.h.j > 0
        ind.b.j <- w.b.j > 0
        N.h.j <- sum(ind.h.j)
        N.b.j <- sum(ind.b.j)
        ind.j <- ind.b.j
        if (h[j] > b[j]) {
          ind.j <- ind.h.j
        }
        ind <- ind.i == "TRUE" | ind.j == "TRUE"
        eN <- sum(ind)
        eY <- y[ind]
        eX <- x[ind]
        W.h.i <- w.h.i[ind]
        W.b.i <- w.b.i[ind]
        W.h.j <- w.h.j[ind]
        W.b.j <- w.b.j[ind]
        eC <- NULL
        if (!is.null(cluster)) {
          eC <- cluster[ind]
        }
        edups <- edupsid <- 0
        if (vce == "nn") {
          edups <- dups[ind]
          edupsid <- dupsid[ind]
        }
        u.i <- (eX - eval[i]) / h[i]
        R.q.i <- matrix(NA, eN, (q + 1))
        for (k in 1:(q + 1)) {
          R.q.i[, k] <- (eX - eval[i])^(k -
            1)
        }
        R.p.i <- R.q.i[, 1:(p + 1)]
        u.j <- (eX - eval[j]) / h[j]
        R.q.j <- matrix(NA, eN, (q + 1))
        for (k in 1:(q + 1)) {
          R.q.j[, k] <- (eX - eval[j])^(k -
            1)
        }
        R.p.j <- R.q.j[, 1:(p + 1)]
        e.p1 <- matrix(0, (q + 1), 1)
        e.p1[p + 2] <- 1
        e.v <- matrix(0, (p + 1), 1)
        e.v[deriv + 1] <- 1
        L.i <- crossprod(R.p.i * W.h.i, u.i^(p + 1))
        invG.q.i <- qrXXinv((sqrt(W.b.i) * R.q.i))
        invG.p.i <- qrXXinv((sqrt(W.h.i) * R.p.i))
        Q.q.i <- t(t(R.p.i * W.h.i) - h[i]^(p + 1) *
          (L.i %*% t(e.p1)) %*% t(t(invG.q.i %*% t(R.q.i)) *
            W.b.i))
        beta.p.i <- invG.p.i %*% crossprod(
          R.p.i * W.h.i,
          eY
        )
        beta.q.i <- invG.q.i %*% crossprod(
          R.q.i * W.b.i,
          eY
        )
        beta.bc.i <- invG.p.i %*% crossprod(Q.q.i, eY)
        L.j <- crossprod(R.p.j * W.h.j, u.j^(p + 1))
        invG.q.j <- qrXXinv((sqrt(W.b.j) * R.q.j))
        invG.p.j <- qrXXinv((sqrt(W.h.j) * R.p.j))
        Q.q.j <- t(t(R.p.j * W.h.j) - h[j]^(p + 1) *
          (L.j %*% t(e.p1)) %*% t(t(invG.q.j %*% t(R.q.j)) *
            W.b.j))
        beta.p.j <- invG.p.j %*% crossprod(
          R.p.j * W.h.j,
          eY
        )
        beta.q.j <- invG.q.j %*% crossprod(
          R.q.j * W.b.j,
          eY
        )
        beta.bc.j <- invG.p.j %*% crossprod(Q.q.j, eY)
        hii.i <- predicts.p.i <- predicts.q.i <- hii.j <- predicts.p.j <- predicts.q.j <- 0
        if (vce == "hc0" | vce == "hc1" | vce == "hc2" |
          vce == "hc3") {
          predicts.p.i <- R.p.i %*% beta.p.i
          predicts.q.i <- R.q.i %*% beta.q.i
          predicts.p.j <- R.p.j %*% beta.p.j
          predicts.q.j <- R.q.j %*% beta.q.j
          if (vce == "hc2" | vce == "hc3") {
            hii.i <- hii.j <- matrix(NA, eN, 1)
            for (k in 1:eN) {
              hii.i[k] <- R.p.i[k, ] %*% invG.p.i %*%
                (R.p.i * W.h.i)[k, ]
              hii.j[k] <- R.p.j[k, ] %*% invG.p.j %*%
                (R.p.j * W.h.j)[k, ]
            }
          }
        }
        res.h.i <- lprobust.res(
          eX, eY, predicts.p.i,
          hii.i, vce, nnmatch, edups, edupsid, p + 1
        )
        if (vce == "nn") {
          res.b.i <- res.h.i
        } else {
          res.b.i <- lprobust.res(
            eX, eY, predicts.q.i,
            hii.i, vce, nnmatch, edups, edupsid, q + 1
          )
        }
        res.h.j <- lprobust.res(
          eX, eY, predicts.p.j,
          hii.j, vce, nnmatch, edups, edupsid, p + 1
        )
        if (vce == "nn") {
          res.b.j <- res.h.j
        } else {
          res.b.j <- lprobust.res(
            eX, eY, predicts.q.j,
            hii.j, vce, nnmatch, edups, edupsid, q + 1
          )
        }
        V.us.i <- factorial(deriv)^2 * invG.p.i %*% t(c(res.h.i) *
          R.p.i * W.h.i)
        V.us.j <- factorial(deriv)^2 * invG.p.j %*% t(c(res.h.j) *
          R.p.j * W.h.j)
        V.rb.i <- factorial(deriv)^2 * invG.p.i %*% t(c(res.b.i) *
          Q.q.i)
        V.rb.j <- factorial(deriv)^2 * invG.p.j %*% t(c(res.b.j) *
          Q.q.j)
        cov.us[i, j] <- (V.us.i %*% t(V.us.j))[deriv +
          1, deriv + 1]
        cov.rb[i, j] <- (V.rb.i %*% t(V.rb.j))[deriv +
          1, deriv + 1]
        cov.us[j, i] <- cov.us[i, j]
        cov.rb[j, i] <- cov.rb[i, j]
      }
    }
  }
  out <- list(Estimate = Estimate, opt = list(
    p = p, q = q,
    deriv = deriv, kernel = kernel.type, n = N, neval = neval,
    bwselect = bwselect
  ), cov.us = cov.us, cov.rb = cov.rb)
  out$call <- match.call()
  class(out) <- "lprobust"
  return(out)
}

# Update nprobust namespace to use custom function
environment(lprobust_cust) <- asNamespace("nprobust")
assignInNamespace("lprobust", lprobust_cust, ns = "nprobust")

get_cv <- function(estimate, y, bc = TRUE, cv = "gcv") {
  #' Convenience function to obtain the cross-validation criterion from the 
  #' estimate matrix returned by lprobust_cust
  n <- nrow(estimate)
  if (bc) {
    y_hat <- estimate[, colnames(estimate) == "tau.bc"]
    A <- estimate[, colnames(estimate) == "A.bc"]
  } else {
    y_hat <- estimate[, colnames(estimate) == "tau.us"]
    A <- estimate[, colnames(estimate) == "A.us"]
  }
  if (cv == "gcv") {
    return(n * sum((y - y_hat)^2) / ((n - sum(A))^2))
  } else if (cv == "ocv") {
    return(mean(((y - y_hat)^2) / ((1 - A)^2)))
  }
}

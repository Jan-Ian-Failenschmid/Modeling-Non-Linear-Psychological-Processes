#' ----------------------------------------------------------------------------#
#' Title: Polynomial regression methods                                        #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 03-09-2024                                                    #
#' -----                                                                       #
#' Last Modified: 13-02-2025                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#

# Polynomial regression with orthogonal polynomials
setClass(
  "method_poly_orth",
  contains = "method"
)

setMethod("fit", "method_poly_orth", function(method, data) {
  model_list <- list()
  p <- 1
  fit <- TRUE
  while (fit) {
    poly_fit <- tryCatch(
      lm(y_obs ~ poly(time, degree = p), data = data),
      error = function(cond) {
        FALSE
      }
    )

    if (class(poly_fit) == "lm") {
      model_list[[p]] <- poly_fit
      p <- p + 1
    } else {
      fit <- FALSE
    }
  }

  gcv_list <- sapply(model_list, function(fit) {
    X <- model.matrix(fit)
    A <- X %*% solve(t(X) %*% X) %*% t(X)
    n <- nrow(X)
    (n * sum((data$y_obs - A %*% data$y_obs)^2)) /
      (n - sum(diag(A)))^2
  })

  fit <- model_list[[which.min(gcv_list)]]

  # Method generics schould always return the adjusted method object
  slot(method, "converged") <- TRUE

  if (slot(method, "converged")) {
    inference <- predict(fit, se.fit = TRUE)

    # Get mean inference
    slot(method, "estimate") <- as.vector(inference$fit)
    slot(method, "ci") <- list(
      ub = as.vector(inference$fit + (qnorm(0.975) * inference$se.fit)),
      lb = as.vector(inference$fit - (qnorm(0.975) * inference$se.fit))
    )

    # Calculate mse
    slot(method, "mse") <- calc_mse(method, data)

    # Calculate gcv
    slot(method, "gcv") <- min(gcv_list)

    # Calculate confidence interval coverage
    slot(method, "ci_coverage") <- ci_test(method, data)
  }

  # Method generics schould always return the adjusted method object
  return(method)
})

#' ----------------------------------------------------------------------------#
#' Title: Simple linear regression methd                                      #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 03-09-2024                                                    #
#' -----                                                                       #
#' Last Modified: 03-09-2024                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#

# SLR
setClass(
  "method_simple",
  contains = "method"
)

setMethod("fit", "method_simple", function(method, data) {
  fit <- lm(y_obs ~ time, data = data)

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
    X <- model.matrix(fit)
    A <- X %*% solve(t(X) %*% X) %*% t(X)
    n <- nrow(X)
    slot(method, "gcv") <-
      (n * sum((data$y_obs - as.vector(inference$fit))^2)) /
        (n - sum(diag(A)))^2

    # Calculate confidence interval coverage
    slot(method, "ci_coverage") <- ci_test(method, data)
  }

  # Method generics schould always return the adjusted method object
  return(method)
})

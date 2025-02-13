#' ----------------------------------------------------------------------------#
#' Title: GAM Methods                                                          #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 11-04-2024                                                    #
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

#' Class and methods for the generalized additive models class
require(mgcv)

# Gams
setClass(
  "method_gam",
  contains = "method"
)

setMethod("fit", "method_gam", function(method, data) {
  # fit <- gam(y_obs ~ s(time, bs = "tp", k = nrow(data)), data = data)
  n <- nrow(data)
  fit <- gam(y_obs ~ s(time, bs = "tp", k = n - 1),
    data = data, method = "ML"
  )


  # Method generics schould always return the adjusted method object
  slot(method, "converged") <- fit$converged

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
    # slot(method, "gcv") <- fit$gcv.ubre
    slot(method, "gcv") <-
      (n * sum((data$y_obs - as.vector(inference$fit))^2)) /
        (n - sum(influence(fit)))^2

    # Calculate confidence interval coverage
    slot(method, "ci_coverage") <- ci_test(method, data)

    # Extract wigglyness parameter
    slot(method, "wiggliness") <- fit$sp
  }

  # Method generics schould always return the adjusted method object
  return(method)
})

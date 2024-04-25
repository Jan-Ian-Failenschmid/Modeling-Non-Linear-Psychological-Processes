#' ----------------------------------------------------------------------------#
#' Title: GAM Methods                                                          #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 11-04-2024                                                    #
#' -----                                                                       #
#' Last Modified: 25-04-2024                                                   #
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
  # ToDo: Have K depend on N for model identifyability.
  fit <- gam(y_obs ~ s(time, bs = "tp", k = nrow(data)), data = data)

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
    slot(method, "gcv") <- fit$gcv.ubre

    # Calculate confidence interval coverage
    slot(method, "ci_coverage") <- ci_test(method, data)
  }

  # Method generics schould always return the adjusted method object
  return(method)
})

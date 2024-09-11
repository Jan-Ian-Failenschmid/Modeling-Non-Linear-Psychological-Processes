#' ----------------------------------------------------------------------------#
#' Title:                                                                      #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 11-09-2024                                                    #
#' -----                                                                       #
#' Last Modified: 11-09-2024                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#

# Polynomial regression
setClass(
  "method_poly",
  contains = "method"
)

setMethod("fit", "method_poly", function(method, data) {
  model_list <- list()
  p <- 1
  fit <- TRUE
  while (fit) {
    poly_fit <- lm(y_obs ~ poly(time, degree = p, raw = TRUE), data = data)

    if (!any(is.na(coef(poly_fit)))) {
      model_list[[p]] <- poly_fit
      p <- p + 1
    } else {
      fit <- FALSE
    }
  }

  gcv_list <- sapply(model_list, function(fit) {
    n <- length(fitted(fit))
    (n * sum((data$y_obs - fitted(fit))^2)) /
      (n - sum(hatvalues(fit)))^2
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

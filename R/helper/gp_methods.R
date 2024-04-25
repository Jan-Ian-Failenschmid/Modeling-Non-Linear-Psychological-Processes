#' ----------------------------------------------------------------------------#
#' Title: GP Methods                                                           #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 12-04-2024                                                    #
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

#' Class and methods for the gaussian process class
require(cmdstanr)

# Gaussian processes
setClass(
  "method_gp",
  contains = "method"
)

setMethod("fit", "method_gp", function(method, data) {
  # Prepare data for STAN
  data <- list(
    N_obs = length(data$time), x_obs = data$time, y_obs = data$y_obs
  )

  # Get stan model
  mod <- quiet(cmdstan_model("./R/stan_files/gp.stan"))

  # Sample from model
  gp_fit <- suppressWarnings(mod$sample(
    data = data,
    seed = 5838298,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    show_messages = FALSE,
    show_exceptions = FALSE,
    adapt_delta = 0.5
  ))

  # Obtain summary measures for the parameters of the gp
  gp_sum <- gp_fit$summary(variables = c("rho", "alpha", "sigma"))[, c(8:10)]

  # Check convergence
  conv <- (
    all(gp_sum$rhat < 1.1) &&
      all(gp_sum$ess_bulk > 500) &&
      all(gp_sum$ess_tail > 500)
  )

  # If not converged rerun with twice the sample size
  if (!conv) {
    gp_fit <- suppressWarnings(mod$sample(
      data = data,
      seed = 5838298,
      chains = 4,
      parallel_chains = 4,
      refresh = 500,
      show_messages = FALSE,
      show_exceptions = FALSE,
      iter_warmup = 1000,
      iter_sampling = 4000,
      adapt_delta = 0.5
    ))

    gp_sum <- gp_fit$summary(variables = c("rho", "alpha", "sigma"))[, c(8:10)]

    conv <- (
      all(gp_sum$rhat < 1.1) &&
        all(gp_sum$ess_bulk > 500) &&
        all(gp_sum$ess_tail > 500)
    )
  }

  slot(method, "converged") <- conv

  if (slot(method, "converged")) {
    # Obtain posterior sample for the expectation of the gaussian process
    posterior_draws <- quiet(
      as.data.frame(gp_fit$draws(
        "f_predict",
        format = "draws_df"
      ))
    )
    posterior_draws <- posterior_draws[
      ,
      seq(1, ncol(posterior_draws) - 3)
    ]

    # Obtain posterior sample for the gaussian process
    posterior_predictive_draws <- quiet(
      as.data.frame(gp_fit$draws(
        "f_post_predict",
        format = "draws_df"
      ))
    )
    posterior_predictive_draws <- posterior_predictive_draws[
      ,
      seq(1, ncol(posterior_predictive_draws) - 3)
    ]

    # Get mean inference
    slot(method, "estimate") <- sapply(posterior_draws, mean)
    slot(method, "ci") <- list(
      ub = sapply(posterior_predictive_draws, quantile, probs = 0.975),
      lb = sapply(posterior_predictive_draws, quantile, probs = 0.025)
    )

    # Calculate mse
    slot(method, "mse") <- calc_mse(method, data)

    # Calculate gcv
    slot(method, "gcv") <- gp_fit$summary(
      variables = "gcv_val",
      mean = "mean"
    )$mean

    # Calculate confidence interval coverage
    slot(method, "ci_coverage") <- ci_test(method, data)
  }

  # Method generics schould always return the adjusted method object
  return(method)
})

#' ----------------------------------------------------------------------------#
#' Title: GP Methods                                                           #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 15-04-2024                                                    #
#' -----                                                                       #
#' Last Modified: 27-05-2024                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#

#' Class and methods for the dynamic modelling analysis method
require(dynr)

# Dynamic modelling
setClass(
  "method_dynm",
  contains = "method"
)

setMethod("fit", "method_dynm", function(method, data) {
  # Fit appropriate dynr model
  iter <- 1
  conv <- FALSE
  delta <- c(0, 0.1, 0.2, 0.3, 0.4)

  # Fit model three times with different starting values in case
  # of non-convergence
  while (iter < 6 && !conv) {
    # Fit the correct dynamic model based on the data generating model
    dynr_fit <- quiet(fit_correct_dynr_model(method, data, delta[iter]))
    # Determine convergence
    conv <- class(dynr_fit) == "dynrCook" &&
      !any(is.na(summary(dynr_fit)$Coefficients)) &&
      dynr_fit$exitflag > 0

    iter <- iter + 1
  }

  slot(method, "converged") <- conv

  if (slot(method, "converged")) {
    # Get mean inference
    y_hat <- dynr_fit@eta_smooth_final[1, ]
    # Get standard error
    se <- sqrt(dynr_fit@error_cov_smooth_final[1, 1, ])

    slot(method, "estimate") <- y_hat
    slot(method, "ci") <- list(
      ub = y_hat + (qnorm(0.975) * se),
      lb = y_hat + (qnorm(0.025) * se)
    )

    # Calculate mse
    slot(method, "mse") <- calc_mse(method, data)

    # Obtain the trace of A
    trace_a <- sum(
      as.vector(dynr_fit@error_cov_smooth_final[1, 1, ]) /
        coef(dynr_fit)["meas_er"]
    )
    # Extract errors
    errors <- data$y_obs - y_hat
    # Calculate gcv
    slot(method, "gcv") <- nrow(data) * sum(errors^2) /
      (nrow(data) - trace_a)^2 # GCV

    # Calculate confidence interval coverage
    slot(method, "ci_coverage") <- ci_test(method, data)
  }

  # Method generics schould always return the adjusted method object
  return(method)
})

fit_correct_dynr_model <- function(method, data, delta = 0) {
  #' Convenience function to fit the appropriate dynamic model with dynr
  # Prepare data
  dynr_dat <- dynr.data(
    data.frame(
      id = 1,
      time = data$time,
      y_obs = data$y_obs
    ),
    id = "id", time = "time", observed = "y_obs"
  )

  if (method@gen_model == "exp_growth") {
    # Measurement model fixed to 1
    meas <- prep.measurement(
      values.load = matrix(1, ncol = 1),
      params.load = matrix("fixed", ncol = 1),
      state.names = c("y"),
      obs.names = c("y_obs")
    )

    # Noise starting values
    mdcov <- prep.noise(
      # Use observed variance as upper bound
      values.latent = diag(var(data$y_obs, na.rm = TRUE), 1),
      params.latent = diag("dyn_er", 1),
      values.observed = diag(var(data$y_obs, na.rm = TRUE), 1),
      params.observed = diag("meas_er", 1)
    )

    # Model equation
    fml_list <- list(
      y ~ yr * ya - yr * y
    )

    slope <- (na.omit(data$y_obs)[length(na.omit(data$y_obs))] -
      na.omit(data$y_obs)[1]) / length(na.omit(data$y_obs))
    max_obs <- max(data$y_obs, na.rm = TRUE)

    # Prepare formula dynamics
    dynm <- prep.formulaDynamics(
      formula = fml_list, startval = c(
        # Start growth rate using the slope between the first and the last value
        yr = c(slope, 1e-6)[which.max(c(slope, 1e-6))] +
          abs(rnorm(1, 0, delta)),
        # Use max observation as starting value
        ya = c(max_obs, 1e-6)[which.max(c(max_obs, 1e-6))] +
          abs(rnorm(1, 0, delta))
      ),
      isContinuousTime = TRUE
    )

    # Prepare initial values
    initial <- prep.initial(
      # Use first observation of y_obs as starting value
      values.inistate = na.omit(data$y_obs)[[1]] + rnorm(1, 0, delta),
      params.inistate = "y_nod",
      # Values fix starting variance to the variance of y_obs as an upper bound
      values.inicov = diag(var(data$y_obs, na.rm = TRUE), 1),
      params.inicov = diag("fixed", 1)
    )

    # Dynr model
    model <- dynr.model(
      dynamics = dynm, measurement = meas, noise = mdcov,
      initial = initial, data = dynr_dat
    )

    # Constraint parameter space to realistic values
    model$lb[c("yr", "ya", "dyn_er", "meas_er", "y_nod")] <- c(
      1e-6, 1e-6, 1e-6, 1e-6, -1e+3
    )

    model$ub[c("yr", "ya", "dyn_er", "meas_er", "y_nod")] <- c(
      1e+3, 1e+3, 1e+3, 1e+3, 1e+3
    )
  } else if (method@gen_model == "log_growth") {
    meas <- prep.measurement(
      values.load = matrix(1, ncol = 1),
      params.load = matrix("fixed", ncol = 1),
      state.names = c("y"),
      obs.names = c("y_obs")
    )

    mdcov <- prep.noise(
      values.latent = diag(var(data$y_obs, na.rm = TRUE), 1),
      params.latent = diag("dyn_er", 1),
      values.observed = diag(var(data$y_obs, na.rm = TRUE), 1),
      params.observed = diag("meas_er", 1)
    )

    fml_list <- list(
      y ~ r * y * (1 - (y / k))
    )

    slope <- (na.omit(data$y_obs)[length(na.omit(data$y_obs))] -
      na.omit(data$y_obs)[1]) / length(na.omit(data$y_obs))
    max_obs <- max(data$y_obs, na.rm = TRUE)

    dynm <- prep.formulaDynamics(
      formula = fml_list, startval = c(
        # Start growth rate using the slope between the first and the last value
        r = c(slope, 1e-6)[which.max(c(slope, 1e-6))] + abs(rnorm(1, 0, delta)),
        # Use max observation as starting value
        k = c(max_obs, 1e-6)[which.max(c(max_obs, 1e-6))] +
          abs(rnorm(1, 0, delta))
      ),
      isContinuousTime = TRUE
    )

    initial <- prep.initial(
      # Prepare initial values
      values.inistate = na.omit(data$y_obs)[[1]] + rnorm(1, 0, delta),
      params.inistate = "y_nod",
      values.inicov = diag(var(data$y_obs, na.rm = TRUE), 1),
      params.inicov = diag("fixed", 1)
    )

    # Dynr model
    model <- dynr.model(
      dynamics = dynm, measurement = meas, noise = mdcov,
      initial = initial, data = dynr_dat
    )

    model$lb[c("r", "k", "dyn_er", "meas_er", "y_nod")] <- c(
      1e-6, 1e-6, 1e-6, 1e-6, -1e+3
    )

    model$ub[c("r", "k", "dyn_er", "meas_er", "y_nod")] <- c(
      1e+3, 1e+3, 1e+3, 1e+3, 1e+3
    )
  } else if (method@gen_model == "damped_oscillator") {
    meas <- prep.measurement(
      values.load = matrix(c(1, 0), 1, 2),
      params.load = matrix(c("fixed", "fixed"), 1, 2),
      state.names = c("y", "v"),
      obs.names = c("y_obs")
    )

    ecov <- prep.noise(
      values.latent = diag(c(var(data$y_obs, na.rm = TRUE), 0), 2),
      params.latent = diag(c("dyn_er", "fixed"), 2),
      values.observed = diag(var(data$y_obs, na.rm = TRUE), 1),
      params.observed = diag("meas_er", 1)
    )

    dynamics <- prep.matrixDynamics(
      values.dyn = matrix(c(
        0, -0.1 - abs(rnorm(1, 0, delta)),
        1, -0.2 - abs(rnorm(1, 0, delta))
      ), 2, 2),
      params.dyn = matrix(c("fixed", "k", "fixed", "c"), 2, 2),
      isContinuousTime = TRUE
    )

    initial <- prep.initial(
      values.inistate = c(
        na.omit(data$y_obs)[[1]] + rnorm(1, 0, delta),
        0 + abs(rnorm(1, 0, delta))
      ),
      params.inistate = c("y_nod", "v_nod"),
      values.inicov = diag(c(var(data$y_obs, na.rm = TRUE), 1), 2),
      params.inicov = diag("fixed", 2)
    )

    model <- dynr.model(
      dynamics = dynamics, measurement = meas, noise = ecov,
      initial = initial, data = dynr_dat
    )
  } else if (method@gen_model == "cusp_catastrophe") {
    meas <- prep.measurement(
      values.load = matrix(c(1, 0, 0), 1, 3),
      params.load = matrix(c("fixed", "fixed", "fixed"), 1, 3),
      state.names = c("y", "b", "v"),
      obs.names = c("y_obs")
    )

    mdcov <- prep.noise(
      values.latent = diag(c(var(data$y_obs, na.rm = TRUE), 0, 0), 3),
      params.latent = diag(c("dyn_er", "fixed", "fixed"), 3),
      values.observed = diag(var(data$y_obs, na.rm = TRUE), 1),
      params.observed = diag("meas_er", 1)
    )

    fml_list <- list(
      y ~ -(4 * y^3 + 2 * a * y + b),
      b ~ v,
      v ~ omega * b
    )

    dynm <- prep.formulaDynamics(
      formula = fml_list, startval = c(
        omega = -0.1 - abs(rnorm(1, 0, delta)),
        a = 0 + rnorm(1, 0, delta)
      ),
      isContinuousTime = TRUE
    )

    initial <- prep.initial(
      # Prepare initial values
      values.inistate = c(
        na.omit(data$y_obs)[[1]] + rnorm(1, 0, delta),
        0 + rnorm(1, 0, delta),
        0 + rnorm(1, 0, delta)
      ),
      params.inistate = c("y_nod", "b_nod", "v_nod"),
      values.inicov = diag(c(var(data$y_obs, na.rm = TRUE), 1, 1), 3),
      params.inicov = diag("fixed", 3)
    )

    model <- dynr.model(
      dynamics = dynm, measurement = meas, noise = mdcov,
      initial = initial, data = dynr_dat
    )

    model$lb[
      c("a", "omega", "dyn_er", "meas_er", "y_nod", "b_nod", "v_nod")
    ] <- c(
      -1e+3, -1e+3, 1e-6, 1e-6, -1e+3, -1e+3, -1e+3
    )

    model$ub[
      c("a", "omega", "dyn_er", "meas_er", "y_nod", "b_nod", "v_nod")
    ] <- c(
      +1e+3, -1e-6, 1e+3, 1e+3, 1e+3, 1e+3, 1e+3
    )
  }

  dynr_fit <- tryCatch(dynr.cook(model),
    error = function(cond) {
      dynr_fit <- "Optimizer Error"
      return(dynr_fit)
    }
  )

  return(dynr_fit)
}

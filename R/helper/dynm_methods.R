#' ----------------------------------------------------------------------------#
#' Title: GP Methods                                                           #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 15-04-2024                                                    #
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

#' Class and methods for the dynamic modelling analysis method
require(dynr)

# Dynamic modelling
setClass(
  "method_dynm",
  contains = "method"
)

setMethod("fit", "method_dynm", function(method, data) {
  # Fit appropriate dynr model
  dynr_fit <- quiet(fit_dynr(method, data))

  # Determine convergence
  if (class(dynr_fit) == "dynrCook") {
    if (any(is.na(summary(dynr_fit)$Coefficients))) {
      conv <- FALSE
    } else {
      conv <- ifelse(dynr_fit$exitflag > 0, TRUE, FALSE)
    }
  } else {
    conv <- FALSE
  }

  # summary(dynr_fit)

  # Method generics schould always return the adjusted method object
  # print(conv)
  slot(method, "fit") <- dynr_fit
  slot(method, "converged") <- conv

  return(method)
})

setMethod(
  "calculate_performance_measures", "method_ssm",
  function(method, data) {
    if (slot(method, "converged")) {
      y_hat <- slot(method, "fit")@eta_smooth_final[1, ]
      se <- sqrt(slot(method, "fit")@error_cov_smooth_final[1, 1, ])
      # Get mean inference
      slot(method, "estimate") <- y_hat
      slot(method, "ci") <- list(
        ub = y_hat + (qnorm(0.975) * se),
        lb = y_hat + (qnorm(0.025) * se)
      )

      # Calculate mse
      slot(method, "mse") <- calc_mse(method, data)

      # ToDo: Test cross-validation for a couple of data sets and see how this
      # works.
      # Calculate gcv
      cv_errors <- sapply(seq_len(nrow(data)), function(x, data, method) {
        target <- data$y_obs[x] # Save target value
        data$y_obs[x] <- NA
        dynr_fit <- fit(method, data)@fit # Fit model with ssm method
        e_i <- target - dynr_fit@eta_smooth_final[1, x] # Calculate cv error
        return(e_i)
      }, data = data, method = method)

      # errors <- data$y_obs - y_hat # Get residuals

      # # Caluclate trace of the influence matrix
      # trace_a <- sum(1 - errors / cv_errors)
      # # Get sample size
      # n <- nrow(data)

      # slot(method, "gcv") <- n * sum(errors^2) / (n - trace_a)^2 # GCV

      slot(method, "gcv") <- mean(cv_errors^2)

      # Calculate confidence interval coverage
      slot(method, "ci_coverage") <- ci_test(method, data)
    } else {
      method <- set_na(method)
    }

    # Method generics schould always return the adjusted method object
    return(method)
  }
)

fit_dynr <- function(method, data) {
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

  if (method@gen_model == "latent_change") {
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

    # Prepare formula dynamics
    dynm <- prep.formulaDynamics(
      formula = fml_list, startval = c(
        # Start growth rate using the slope between the first and the last value
        yr = (na.omit(data$y_obs)[length(na.omit(data$y_obs))] -
          na.omit(data$y_obs)[1]) / length(na.omit(data$y_obs)),
        # Use max observation as starting value
        ya = max(data$y_obs, na.rm = TRUE)
      ),
      isContinuousTime = TRUE
    )

    # Prepare initial values
    initial <- prep.initial(
      # Use first observation of y_obs as starting value
      values.inistate = na.omit(data$y_obs)[[1]],
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

    dynm <- prep.formulaDynamics(
      formula = fml_list, startval = c(
        r = (na.omit(data$y_obs)[length(na.omit(data$y_obs))] -
          na.omit(data$y_obs)[1]) / length(na.omit(data$y_obs)),
        k = max(data$y_obs, na.rm = TRUE)
      ),
      isContinuousTime = TRUE
    )

    initial <- prep.initial(
      # Prepare initial values
      values.inistate = na.omit(data$y_obs)[[1]],
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
        0, -0.1, 1,
        -0.2
      ), 2, 2),
      params.dyn = matrix(c("fixed", "k", "fixed", "c"), 2, 2),
      isContinuousTime = TRUE
    )

    initial <- prep.initial(
      values.inistate = c(na.omit(data$y_obs)[[1]], 0),
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
      formula = fml_list, startval = c(omega = -0.1, a = 0),
      isContinuousTime = TRUE
    )

    initial <- prep.initial(
      # Prepare initial values
      values.inistate = c(na.omit(data$y_obs)[[1]], 0, 0),
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

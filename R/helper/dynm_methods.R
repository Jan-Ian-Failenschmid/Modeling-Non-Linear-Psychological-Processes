#' ----------------------------------------------------------------------------#
#' Title: GP Methods                                                           #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 15-04-2024                                                    #
#' -----                                                                       #
#' Last Modified: 06-09-2024                                                   #
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
require(future.apply)

# Dynamic modelling
setClass(
  "method_dynm",
  contains = "method"
)

setMethod("fit", "method_dynm", function(method, data) {
  n <- 5 * 4
  if (method@gen_model == "exp_growth") {
    start_val <- data.frame(
      b = runif(n, -1, 1e-6),
      a = runif(n, 1e-6, 1),
      dyn_er = runif(n, 1e-6, var(data$y_obs, na.rm = TRUE)),
      meas_er = runif(n, 1e-6, var(data$y_obs, na.rm = TRUE)),
      y_nod = runif(n, -10, 10)
    )

    dynr_fit <- quiet(fit_dynm(fit_dynm_exp_growth, start_val, data))
  } else if (method@gen_model == "log_growth") {
    start_val <- data.frame(
      r = runif(n, 1e-6, 5),
      k = runif(n, 1e-6, 50),
      dyn_er = runif(n, 1e-6, var(data$y_obs, na.rm = TRUE)),
      meas_er = runif(n, 1e-6, var(data$y_obs, na.rm = TRUE)),
      y_nod = runif(n, 1e-6, 10)
    )
    dynr_fit <- quiet(fit_dynm(fit_dynm_log_growth, start_val, data))
  } else if (method@gen_model == "cusp_catastrophe") {
    start_val <- data.frame(
      omega = runif(n, -5, 1e-6),
      a = runif(n, -10, 10),
      dyn_er = runif(n, 0, var(data$y_obs, na.rm = TRUE)),
      meas_er = runif(n, 0, var(data$y_obs, na.rm = TRUE)),
      y_nod = runif(n, -10, 10),
      b_nod = runif(n, -50, 50),
      v_nod = runif(n, -10, 10)
    )

    dynr_fit <- quiet(fit_dynm(fit_dynm_cusp_catas, start_val, data))
  } else if (method@gen_model == "damped_oscillator") {
    start_val <- data.frame(
      k = runif(n, -1, 1e-6),
      c = runif(n, -1, 1e-6),
      dyn_er = runif(n, 1e-6, var(data$y_obs, na.rm = TRUE)),
      meas_er = runif(n, 1e-6, var(data$y_obs, na.rm = TRUE)),
      y_nod = runif(n, -10, 10),
      v_nod = runif(n, -10, 10)
    )

    dynr_fit <- quiet(fit_dynm(fit_dynm_damp_osc, start_val, data))
  } else {
    stop("Parametric model not defined!")
  }

  # Determine convergence
  conv <- class(dynr_fit) == "dynrCook" &&
    !any(is.na(summary(dynr_fit)$Coefficients)) &&
    dynr_fit$exitflag > 0 &&
    max(dynr_fit@error_cov_smooth_final[1, 1, ]) >= 1e-3

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

fit_dynm <- function(fit_fun, start, data) {
  for (i in 1:2) {
    if (i == 1) {
      fit <- future_apply(start, 1, fit_fun,
        data = data,
        simplify = FALSE,
        future.seed = TRUE
      )
    } else if (i == 2) {
      new_start <- do.call(
        rbind,
        lapply(
          fit[order(dev, decreasing = FALSE)[1:(nrow(start) / 4)]],
          function(fit) {
            if (class(fit) == "dynrCook") {
              coef(fit)
            } else {
              start[sample(seq_len(nrow(start)), 1), ]
            }
          }
        )
      )

      new_start <- apply(
        new_start[rep(seq_len(nrow(new_start)), each = 2), ], 2,
        jitter
      )

      fit <- list(
        fit,
        future_apply(new_start, 1, fit_fun,
          data = data,
          simplify = FALSE,
          future.seed = TRUE
        )
      )
    }

    dev <- sapply(fit, function(x) {
      if (class(x) == "dynrCook") {
        if (any(is.na(x$standard.errors))) {
          1e+9
        } else {
          deviance(x)
        }
      } else {
        1e+9
      }
    })
  }
  return(fit[which.min(dev)][[1]])
}

fit_dynm_log_growth <- function(start, data) {
  dynr_dat <- dynr.data(
    data.frame(
      id = 1,
      time = data$time,
      y_obs = data$y_obs
    ),
    id = "id", time = "time", observed = "y_obs"
  )

  meas <- prep.measurement(
    values.load = matrix(1, ncol = 1),
    params.load = matrix("fixed", ncol = 1),
    state.names = c("y"),
    obs.names = c("y_obs")
  )

  mdcov <- prep.noise(
    values.latent = diag(start[3], 1),
    params.latent = diag("dyn_er", 1),
    values.observed = diag(start[4], 1),
    params.observed = diag("meas_er", 1)
  )

  fml_list <- list(
    y ~ r * y * (1 - (y / k))
  )

  dynm <- prep.formulaDynamics(
    formula = fml_list,
    startval = unlist(start[1:2]),
    isContinuousTime = TRUE
  )

  initial <- prep.initial(
    # Prepare initial values
    values.inistate = start[5],
    params.inistate = "y_nod",
    # diag(var(data$y_obs, na.rm = TRUE), 1),
    values.inicov = diag(100, 1),
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

  dynr_fit <- tryCatch(dynr.cook(model),
    error = function(cond) {
      dynr_fit <- "Optimizer Error"
      return(dynr_fit)
    }
  )
}

fit_dynm_exp_growth <- function(start, data) {
  dynr_dat <- dynr.data(
    data.frame(
      id = 1,
      time = data$time,
      y_obs = data$y_obs
    ),
    id = "id", time = "time", observed = "y_obs"
  )

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
    values.latent = diag(start[3], 1),
    params.latent = diag("dyn_er", 1),
    values.observed = diag(start[4], 1),
    params.observed = diag("meas_er", 1)
  )

  dynm <- prep.matrixDynamics(
    values.dyn = matrix(start[1]),
    params.dyn = matrix("b"),
    values.int = start[2],
    params.int = "a",
    isContinuousTime = TRUE
  )

  # Prepare initial values
  initial <- prep.initial(
    # Use first observation of y_obs as starting value
    values.inistate = start[5],
    params.inistate = "y_nod",
    # Values fix starting variance to the variance of y_obs as an upper bound
    values.inicov = diag(100, 1),
    params.inicov = diag("fixed", 1)
  )

  # Dynr model
  model <- dynr.model(
    dynamics = dynm, measurement = meas, noise = mdcov,
    initial = initial, data = dynr_dat
  )

  # Constraint parameter space to realistic values
  model$lb[c("b", "a", "dyn_er", "meas_er", "y_nod")] <- c(
    -1e+3, 1e-6, 1e-6, 1e-6, -1e+3
  )

  model$ub[c("b", "a", "dyn_er", "meas_er", "y_nod")] <- c(
    -1e-6, 1e+3, 1e+3, 1e+3, 1e+3
  )

  dynr_fit <- tryCatch(dynr.cook(model),
    error = function(cond) {
      dynr_fit <- "Optimizer Error"
      return(dynr_fit)
    }
  )
}

fit_dynm_damp_osc <- function(start, data) {
  dynr_dat <- dynr.data(
    data.frame(
      id = 1,
      time = data$time,
      y_obs = data$y_obs
    ),
    id = "id", time = "time", observed = "y_obs"
  )

  meas <- prep.measurement(
    values.load = matrix(c(1, 0), 1, 2),
    params.load = matrix(c("fixed", "fixed"), 1, 2),
    state.names = c("y", "v"),
    obs.names = c("y_obs")
  )

  ecov <- prep.noise(
    values.latent = diag(c(start[3], 0), 2),
    params.latent = diag(c("dyn_er", "fixed"), 2),
    values.observed = diag(start[4], 1),
    params.observed = diag("meas_er", 1)
  )

  dynamics <- prep.matrixDynamics(
    values.dyn = matrix(c(
      0, start[1],
      1, start[2]
    ), 2, 2),
    params.dyn = matrix(c("fixed", "k", "fixed", "c"), 2, 2),
    isContinuousTime = TRUE
  )

  initial <- prep.initial(
    values.inistate = c(
      start[5],
      start[6]
    ),
    params.inistate = c("y_nod", "v_nod"),
    values.inicov = diag(c(100, 10), 2),
    params.inicov = diag("fixed", 2)
  )

  model <- dynr.model(
    dynamics = dynamics, measurement = meas, noise = ecov,
    initial = initial, data = dynr_dat
  )

  dynr_fit <- tryCatch(dynr.cook(model),
    error = function(cond) {
      dynr_fit <- "Optimizer Error"
      return(dynr_fit)
    }
  )

  if (class(dynr_fit) == "dynrCook") {
    c(deviance(dynr_fit), coef(dynr_fit))
  } else {
    c(1e10, start)
  }
}


fit_dynm_cusp_catas <- function(start, data) {
  dynr_dat <- dynr.data(
    data.frame(
      id = 1,
      time = data$time,
      y_obs = data$y_obs
    ),
    id = "id", time = "time", observed = "y_obs"
  )

  meas <- prep.measurement(
    values.load = matrix(c(1, 0, 0), 1, 3),
    params.load = matrix(c("fixed", "fixed", "fixed"), 1, 3),
    state.names = c("y", "b", "v"),
    obs.names = c("y_obs")
  )

  mdcov <- prep.noise(
    values.latent = diag(c(start[[3]], 0, 0), 3),
    params.latent = diag(c("dyn_er", "fixed", "fixed"), 3),
    values.observed = diag(start[[4]], 1),
    params.observed = diag("meas_er", 1)
  )

  fml_list <- list(
    y ~ -(4 * y^3 + 2 * a * y + b),
    b ~ v,
    v ~ omega * b
  )

  dynm <- prep.formulaDynamics(
    formula = fml_list, startval = c(
      omega = start[[1]],
      a = start[[2]]
    ),
    isContinuousTime = TRUE
  )

  initial <- prep.initial(
    # Prepare initial values
    values.inistate = c(
      start[[5]],
      start[[6]],
      start[[7]]
    ),
    params.inistate = c("y_nod", "b_nod", "v_nod"),
    values.inicov = diag(c(100, 10, 10), 3),
    params.inicov = diag("fixed", 3)
  )

  model <- dynr.model(
    dynamics = dynm, measurement = meas, noise = mdcov,
    initial = initial, data = dynr_dat
  )

  # model$lb[
  #   c("a", "omega", "dyn_er", "meas_er", "y_nod", "b_nod", "v_nod")
  # ] <- c(
  #   -1e+3, -1e+3, 1e-6, 1e-6, -1e+3, -1e+3, -1e+3
  # )

  # model$ub[
  #   c("a", "omega", "dyn_er", "meas_er", "y_nod", "b_nod", "v_nod")
  # ] <- c(
  #   +1e+3, -1e-6, 1e+3, 1e+3, 1e+3, 1e+3, 1e+3
  # )

  dynr_fit <- tryCatch(dynr.cook(model),
    error = function(cond) {
      dynr_fit <- "Optimizer Error"
      return(dynr_fit)
    }
  )
}

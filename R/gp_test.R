#' ----------------------------------------------------------------------------#
#' Title:                                                                      #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 14-06-2024                                                    #
#' -----                                                                       #
#' Last Modified: 26-06-2024                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#
library(data.table)
library(cmdstanr)
library(mgcv)

cov_fun <- function(x, y) {
  cov_fun_sp <- function(x, y) {
    dist <- abs(x - y)
    v <- min(x, y)
    return((dist * v^2) / 2 + (v^3 / 3))
  }
  cov_mat <- matrix(0, nrow = length(x), ncol = length(y))
  for (n1 in seq_len(length(x))) {
    for (n2 in seq_len(length(y))) {
      cov_mat[n1, n2] <- cov_fun_sp(x[n1], y[n2])
    }
  }
  return(cov_mat)
}

mse <- function(pred, val) {
  mean((pred - val)^2)
}

# Load data ------
load("R/data/simulation_data_17_05_2024_06_59.Rdata")
load("R/data/simulation_results_17_05_2024_06_59.Rdata")

par(mfrow = c(3, 2))

for (fig in c(seq_len(6))) {
  if (fig == 1) {
    # Generate data
    data <- data.frame(time = seq(0, 1, length.out = 202)[2:201])
    cov_mat <- cov_fun(data$time, data$time)
    data$y <- 1 + data$time + 5 * chol(cov_mat) %*% rnorm(200, 0, 1)
    data$y_obs <- data$y + rnorm(200, 0, 1)
    main <- "GP"
  } else if (fig == 2) {
    # Generate data
    data <- data.frame(time = seq(0, 1, length.out = 202)[2:201])
    data$y <- 1 + data$time + sin(10 * data$time)
    data$y_obs <- data$y + rnorm(200, 0, 1)
    main <- "Sin + Trend"
  } else if (fig == 3) {
    # Generate data
    data <- sim$dat[[ceiling(14995 / 4)]]
    main <- "K = nrow()"
  } else if (fig == 4) {
    # Generate data
    data <- sim$dat[[ceiling(14995 / 4)]]
    main <- "K = nrow() - 1"
  } else if (fig == 5) {
    # Generate data
    data <- sim$dat[[ceiling(12255 / 4)]]
    main <- "K = nrow()"
  } else if (fig == 6) {
    # Generate data
    data <- sim$dat[[ceiling(12255 / 4)]]
    main <- "K = -1"
  }

  plot(data$time, data$y_obs, main = main)
  lines(data$time, data$y)

  # Fit gp --------
  stan_data <- list(
    N_obs = length(data$time), x_obs = data$time / max(data$time),
    y_obs = as.vector(data$y_obs),
    x_pred = data$time / max(data$time), N_pred = length(data$time)
  )

  mod <- cmdstan_model("./R/stan_files/gp_test.stan")
  gp_fit <- mod$optimize(data = stan_data)
  gp_pred <- gp_fit$summary("f")$estimate
  lines(data$time, gp_pred, col = "blue")

  # Fit old gp -------
  mod1 <- cmdstan_model("./R/stan_files/gp.stan")
  gp_fit1 <- mod1$sample(data = stan_data, parallel_chains = 4)
  gp_pred1 <- gp_fit1$summary("f_predict")$mean
  lines(data$time, gp_pred1, col = "green")

  # Fit gam --------
  if (fig %in% c(1, 2, 3, 5)) {
    form = formula(data$y_obs ~ s(data$time, bs = "tp", k = nrow(data)))
  } else if (fig == 4) {
    form = formula(data$y_obs ~ s(data$time, bs = "tp", k = nrow(data) - 1))
  } else if (fig == 6) {
    form = formula(data$y_obs ~ s(data$time, bs = "tp"))
  }

  gam_fit <- gam(form)
  gam_pred <- fitted(gam_fit)
  lines(data$time, gam_pred, col = "red")
}

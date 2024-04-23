#' ----------------------------------------------------------------------------#
#' Title: Method Illustration Figure                                           #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 04-04-2024                                                    #
#' -----                                                                       #
#' Last Modified: 18-04-2024                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#



### Set-up ---------------------------------------------------------------------
pacman::p_load(nprobust, mgcv, cmdstanr)
# Load functions
invisible(sapply(
  c(paste0("./R/helper/", dir(path = "./R/helper"))),
  source
))

plot_local_polynomial <- function(
    x, y, x_eval, h, p = 1,
    K = dnorm, threshold, alpha = 1) {
  # Create X matrix of the local polinomial smoother (x_i - x)^p
  X <- as.matrix(sapply(0:p, function(p, x, x_eval) (x - x_eval)^p,
    x = x, x_eval = x_eval
  ))

  # Calculate weights for the distances between each x points and x x.eval
  W <- diag(1, length(x)) * (K((x_eval - x) / h))

  # Solve local polynomial equation to get to influence vector
  beta <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y

  # Make predictions for the local polynomial
  pred <- X %*% beta

  lines(x[which(diag(W) > threshold)], pred[which(diag(W) > threshold)],
    col = rgb(red = 1, green = 0, blue = 0, alpha = alpha)
  )

  points(x[x_eval], pred[x_eval],
    pch = 19,
    col = rgb(red = 1, green = 0, blue = 0, alpha = alpha),
    cex = 2
  )

  return(NULL)
}

set.seed(41)
fig_path <- c("./figures/")
cex <- 2.5
# Simulate data ----------------------------------------------------------------
damp_osc <- new("gen_model",
  time = 200,
  model_name = "damped_oscillator",
  model_type = "DE",
  delta = 1 / 30,
  pars = list(b = 0.01, omega = 0.1),
  start = list(
    formula(y ~ 2),
    formula(v ~ 0)
  ),
  state_eq = list(
    formula(y ~ v),
    formula(v ~ -2 * b * v - (omega^2) * y)
  ),
  meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, 1)))
)

dat <- get_tsm_data(sim_tsm(damp_osc), T, T, T)

### Fits and Plots -------------------------------------------------------------
# Local polynomial
lp_fit <- lprobust(dat$y_obs, dat$time,
  eval = dat$time, p = 3,
  kernel = "gau", bwselect = "imse-dpi",
  bwcheck = 0
)

png(
  file = paste0(fig_path, "locpol_demonstration.png"),
  width = 1960, height = 1080
)

plot(dat$time, dat$y_obs,
  main = "Local polynomial regression demonstration",
  xlab = "Time", ylab = "Y",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex
) # Data points
lines(dat$time, dat$y, lty = 2, lwd = 3) # State function
lines(dat$time, lp_fit$Estimate[, 5], lwd = 3) # Estimated line
# Local polynomials
for (eval in c(50, 100, 150)) {
  plot_local_polynomial(
    x = dat$time, y = dat$y_obs, x_eval = eval,
    h = lp_fit$Estimate[1, 2], p = lp_fit$opt$p, threshold = 0.02
  )
}
# Legend
legend(150, 4,
  legend = c("Process", "Prediction", "Local polynomials"),
  col = c("black", "black", "red"), lty = c(2, 1, 1), cex = cex
)

dev.off()

# GAM
gam_fit <- gam(y_obs ~ s(time, bs = "tp", k = 10), data = dat)
predmat <- predict(gam_fit, type = "lpmatrix")
predmat <- as.data.frame(t(t(predmat) * gam_fit$coefficients))

png(
  file = paste0(fig_path, "gam_demonstration.png"),
  width = 1960, height = 1080
)

plot(dat$time, dat$y_obs,
  main = "Generalized additive model demonstration",
  xlab = "Time", ylab = "Y",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex
) # Data points
lines(dat$time, dat$y, lty = 2, lwd = 3) # State functions
lines(dat$time, fitted(gam_fit), lwd = 3) # Predictions
# Basis function values
for (i in seq_len(ncol(predmat))) {
  lines(
    x = dat$time, y = predmat[, i],
    col = rgb(red = 1, green = 0, blue = 0, alpha = 1)
  )
}
# Legend
legend(150, 4,
  legend = c("Process", "Prediction", "Basis function value"),
  col = c("black", "black", "red"), lty = c(2, 1, 1), cex = cex
)

dev.off()

# GP
data <- list(
  N_obs = length(dat$time), x_obs = dat$time, y_obs = dat$y_obs,
  N_predict = length(dat$time),
  x_predict = dat$time
)

mod <- cmdstan_model("./R/stan_files/gp3.stan")

gp_fit <- mod$sample(
  data = data,
  seed = 5838298,
  chains = 4,
  parallel_chains = 4,
  refresh = 500
)

png(
  file = paste0(fig_path, "gp_demonstration.png"),
  width = 1960, height = 1080
)

# Posterior predictive draws
gp_post_pred(gp_fit,
  f_name = "f_post_predict", time = dat$time,
  obs = dat$y_obs, state = rep(-20, 200), alpha = 0.2,
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex,
  main = "Gaussian process demonstration",
  xlab = "Time", ylab = "Y",
)
lines(dat$time, dat$y, lty = 2, lwd = 3)
lines(x = dat$time, y = gp_fit$summary("f_post_predict")$mean, lwd = 3)
legend(150, 4,
  legend = c("Process", "Prediction", "Posterior function draws"),
  col = c("black", "black", "red"), lty = c(2, 1, 1), cex = cex
)

dev.off()

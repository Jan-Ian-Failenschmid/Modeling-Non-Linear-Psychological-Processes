#' ----------------------------------------------------------------------------#
#' Title: Method Illustration Figure                                           #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 04-04-2024                                                    #
#' -----                                                                       #
#' Last Modified: 09-02-2025                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#



### Set-up ---------------------------------------------------------------------
library(nprobust)
library(mgcv)
library(cmdstanr)
library(MASS)
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

  points(x[which(x == x_eval)], pred[which(x == x_eval)],
    pch = 19,
    col = rgb(red = 1, green = 0, blue = 0, alpha = alpha),
    cex = 2
  )

  return(NULL)
}

set.seed(41)
fig_path <- c("./figures/")
cex <- 3.5
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

dat <- get_tsm_data(sim_tsm(damp_osc))

# Fits and Plots ---------------------------------------------------------------
## Local polynomial ------
png(
  file = paste0(fig_path, "locpol_demonstration.png"),
  width = 1960, height = 1080
)

layout_mat <- matrix(
  c(
    2, 2, 3, 3, 4, 4,
    5, 5, 6, 6, 7, 7,
    8, 8, 9, 9, 10, 10,
    1, 1, 1, 1, 1, 1
  ), 4, 6,
  byrow = TRUE
)
layout(mat = layout_mat, heights = c(0.3, 0.3, 0.3, 0.1))

par(lwd = 2.5)

# Top right plot
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
# Legend
legend(
  x = "bottom",
  legend = c("Process", "Prediction", "Local polynomials"),
  col = c("black", "black", "red"), lty = c(2, 1, 1), cex = cex,
  lwd = rep(2.5, 3),
  xpd = TRUE, horiz = TRUE, seg.len = 1, bty = "n"
)

lp_fit <- lprobust(dat$y_obs, dat$time,
  eval = dat$time[!(dat$time %% 2)], p = 3,
  kernel = "epa", bwselect = "imse-dpi",
  bwcheck = 0, imsegrid = 200
)

# Top left plot
plot(dat$time, dat$y_obs,
  main = "a) Fitted local polynomial regression",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE
) # Data points

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

lines(dat$time, dat$y, lty = 2, lwd = 2.5) # State function
lines(dat$time[!(dat$time %% 2)], lp_fit$Estimate[, 5], lwd = 2.5)


epanechnikov_kernel <- function(x) {
  ifelse(abs(x) <= 1, 0.75 * (1 - x^2), 0)
}

# Estimated line
# Local polynomials
for (eval in c(98, 100, 102)) {
  plot_local_polynomial(
    x = dat$time, y = dat$y_obs, x_eval = eval,
    h = lp_fit$Estimate[1, 2], p = lp_fit$opt$p, threshold = 0.02,
    K = epanechnikov_kernel
  )
}

# Rectangle
rect(
  xleft = 95, xright = 105,
  ybottom = lp_fit$Estimate[, 5][which(lp_fit$Estimate[, 1] == 100)] - .2,
  ytop = lp_fit$Estimate[, 5][which(lp_fit$Estimate[, 1] == 100)] + .2,
  lwd = 4
)

# Top mid plot
plot(dat$time, dat$y_obs,
  main = "b) Zoomed in local polynomial regression",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE,
  xlim = c(95, 105),
  ylim = c(lp_fit$Estimate[, 5][which(lp_fit$Estimate[, 1] == 100)] - .2, lp_fit$Estimate[, 5][which(lp_fit$Estimate[, 1] == 100)] + .2)
) # Data points

axis(2, at = seq(-10, 10, .2), cex.axis = cex)
axis(1,
  at = c(95, 97.5, 100, 12.5, 105),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

lines(dat$time, dat$y, lty = 2, lwd = 2.5) # State function
lines(dat$time[!(dat$time %% 2)], lp_fit$Estimate[, 5], lwd = 2.5)
# Estimated line
# Local polynomials
for (eval in c(98, 100, 102)) {
  plot_local_polynomial(
    x = dat$time, y = dat$y_obs, x_eval = eval,
    h = lp_fit$Estimate[1, 2], p = lp_fit$opt$p, threshold = 0.02,
    K = epanechnikov_kernel
  )
}

# Top right plot
loc <- seq(-1.5, 1.5, 0.01)
weight <- sapply((loc - 0), epanechnikov_kernel)

plot(loc, weight,
  type = "l",
  main = "c) Epanechnikov Kernel",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE
) # Data points

axis(2, at = seq(-10, 10, 0.5), cex.axis = cex)
axis(1,
  at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
  labels = c("", "-1", "", "", "", "1", ""),
  cex.axis = cex, padj = 1
)
title(xlab = expression(paste(Delta, "t / h")), line = 4, cex.lab = cex)

# Mid left plot
lp_fit <- lprobust(dat$y_obs, dat$time,
  eval = dat$time[!(dat$time %% 2)], p = 3,
  kernel = "epa", h = 5,
  bwcheck = 0, imsegrid = 200
)

plot(dat$time, dat$y_obs,
  main = "d) Bandwidth is too small",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE
) # Data points

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

lines(dat$time, dat$y, lty = 2, lwd = 2.5) # State function
lines(dat$time[!(dat$time %% 2)], lp_fit$Estimate[, 5], lwd = 2.5)

# Estimated line
# Local polynomials
for (eval in c(50, 100, 150)) {
  plot_local_polynomial(
    x = dat$time, y = dat$y_obs, x_eval = eval,
    h = lp_fit$Estimate[1, 2], p = lp_fit$opt$p, threshold = 0.02,
    K = epanechnikov_kernel
  )
}


# Mid center plot
lp_fit <- lprobust(dat$y_obs, dat$time,
  eval = dat$time[!(dat$time %% 2)], p = 3,
  kernel = "epa", bwselect = "imse-dpi",
  bwcheck = 0, imsegrid = 200
)

plot(dat$time, dat$y_obs,
  main = "e) Bandwidth is optimized",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE
) # Data points

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

lines(dat$time, dat$y, lty = 2, lwd = 2.5) # State function
lines(dat$time[!(dat$time %% 2)], lp_fit$Estimate[, 5], lwd = 2.5)
# Estimated line
# Local polynomials
for (eval in c(50, 100, 150)) {
  plot_local_polynomial(
    x = dat$time, y = dat$y_obs, x_eval = eval,
    h = lp_fit$Estimate[1, 2], p = lp_fit$opt$p, threshold = 0.02,
    K = epanechnikov_kernel
  )
}


# Mid right plot
lp_fit <- lprobust(dat$y_obs, dat$time,
  eval = dat$time[!(dat$time %% 2)], p = 3,
  kernel = "epa", h = 75,
  bwcheck = 0, imsegrid = 200
)

plot(dat$time, dat$y_obs,
  main = "f) Bandwidth is too large",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE
) # Data points

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

lines(dat$time, dat$y, lty = 2, lwd = 2.5) # State function
lines(dat$time[!(dat$time %% 2)], lp_fit$Estimate[, 5], lwd = 2.5)
# Estimated line
# Local polynomials
for (eval in c(50, 100, 150)) {
  plot_local_polynomial(
    x = dat$time, y = dat$y_obs, x_eval = eval,
    h = lp_fit$Estimate[1, 2], p = lp_fit$opt$p, threshold = 0.02,
    K = epanechnikov_kernel
  )
}

# Bottom left plot
lp_fit <- lprobust(dat$y_obs, dat$time,
  eval = dat$time[!(dat$time %% 2)], p = 1,
  kernel = "epa", bwselect = "imse-dpi",
  bwcheck = 0, imsegrid = 200
)

plot(dat$time, dat$y_obs,
  main = "g) Local linear",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE
) # Data points

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

lines(dat$time, dat$y, lty = 2, lwd = 2.5) # State function
lines(dat$time[!(dat$time %% 2)], lp_fit$Estimate[, 5], lwd = 2.5)
# Estimated line
# Local polynomials
for (eval in c(50, 100, 150)) {
  plot_local_polynomial(
    x = dat$time, y = dat$y_obs, x_eval = eval,
    h = lp_fit$Estimate[1, 2], p = lp_fit$opt$p, threshold = 0.02,
    K = epanechnikov_kernel
  )
}

# Mid center plot
lp_fit <- lprobust(dat$y_obs, dat$time,
  eval = dat$time[!(dat$time %% 2)], p = 2,
  kernel = "epa", bwselect = "imse-dpi",
  bwcheck = 0, imsegrid = 200
)

plot(dat$time, dat$y_obs,
  main = "h) Local quadratic",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE
) # Data points

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

lines(dat$time, dat$y, lty = 2, lwd = 2.5) # State function
lines(dat$time[!(dat$time %% 2)], lp_fit$Estimate[, 5], lwd = 2.5)
# Estimated line
# Local polynomials
for (eval in c(50, 100, 150)) {
  plot_local_polynomial(
    x = dat$time, y = dat$y_obs, x_eval = eval,
    h = lp_fit$Estimate[1, 2], p = lp_fit$opt$p, threshold = 0.02,
    K = epanechnikov_kernel
  )
}

# Mid right plot
lp_fit <- lprobust(dat$y_obs, dat$time,
  eval = dat$time[!(dat$time %% 2)], p = 4,
  kernel = "epa", bwselect = "imse-dpi",
  bwcheck = 0, imsegrid = 200
)

plot(dat$time, dat$y_obs,
  main = "i) Local quartic",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE
) # Data points

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

lines(dat$time, dat$y, lty = 2, lwd = 2.5) # State function
lines(dat$time[!(dat$time %% 2)], lp_fit$Estimate[, 5], lwd = 2.5)
# Estimated line
# Local polynomials
for (eval in c(50, 100, 150)) {
  plot_local_polynomial(
    x = dat$time, y = dat$y_obs, x_eval = eval,
    h = lp_fit$Estimate[1, 2], p = lp_fit$opt$p, threshold = 0.02,
    K = epanechnikov_kernel
  )
}

dev.off()

## GAM -----
png(
  file = paste0(fig_path, "gam_demonstration.png"),
  width = 1960, height = 1080
)

layout_mat <- matrix(
  c(
    2, 2, 3, 3,
    4, 4, 5, 5,
    6, 6, 7, 7,
    1, 1, 1, 1
  ), 4, 4,
  byrow = TRUE
)
layout(mat = layout_mat, heights = c(0.3, 0.3, 0.3, 0.1))

par(lwd = 2.5)

# Legend
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
# Legend
legend(
  x = "bottom",
  legend = c("Process", "Prediction", "Basis function"),
  col = c("black", "black", "red"), lty = c(2, 1, 1), cex = cex,
  lwd = rep(2.5, 3),
  xpd = TRUE, horiz = TRUE, seg.len = 1, bty = "n"
)

par(lwd = 2.5)

# Top left plot
gam_fit <- gam(y_obs ~ s(time, bs = "tp", k = 10), data = dat)
predmat <- predict(gam_fit, type = "lpmatrix")
predmat <- as.data.frame(t(t(predmat) * gam_fit$coefficients))

plot(dat$time, dat$y_obs,
  type = "n",
  main = "a) Fitted GAM with weighted basis functions",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE
) # Data points

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

# Basis function values
for (i in seq_len(ncol(predmat))) {
  lines(
    x = dat$time, y = predmat[, i],
    col = rgb(red = 1, green = 0, blue = 0, alpha = 1)
  )
}
lines(dat$time, dat$y, lty = 2, lwd = 3) # State functions
lines(dat$time, fitted(gam_fit), lwd = 3) # Predictions

points(dat$time, dat$y_obs)

# Top right plot
gam_fit <- gam(y_obs ~ s(time, bs = "tp", k = 10), data = dat)
predmat <- as.data.frame(predict(gam_fit, type = "lpmatrix"))

plot(dat$time, dat$y_obs,
  type = "n",
  main = "b) Unweighted basis functions",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE
) # Data points

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

# Basis function values
for (i in seq_len(ncol(predmat))) {
  lines(
    x = dat$time, y = predmat[, i],
    col = rgb(red = 1, green = 0, blue = 0, alpha = 1)
  )
}
lines(dat$time, dat$y, lty = 2, lwd = 3) # State functions
# lines(dat$time, fitted(gam_fit), lwd = 3) # Predictions

points(dat$time, dat$y_obs)

# Top right plot
# gam_fit <- gam(y_obs ~ s(time, bs = "cr", k = 10), data = dat)
# predmat <- predict(gam_fit, type = "lpmatrix")
# predmat <- as.data.frame(t(t(predmat) * gam_fit$coefficients))

# plot(dat$time, dat$y_obs,
#   type = "n",
#   main = "Fitted GAM with different basis functions",
#   xlab = "", ylab = "",
#   cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE
# ) # Data points

# axis(2, at = seq(-10, 10, 2), cex.axis = cex)
# axis(1,
#   at = c(0, 50, 100, 150, 200),
#   labels = c("", "", "", "", ""),
#   cex.axis = cex, padj = 1
# )
# title(xlab = "Time", line = 4, cex.lab = cex)

# # Basis function values
# for (i in seq_len(ncol(predmat))) {
#   lines(
#     x = dat$time, y = predmat[, i],
#     col = rgb(red = 1, green = 0, blue = 0, alpha = 1)
#   )
# }
# lines(dat$time, dat$y, lty = 2, lwd = 3) # State functions
# lines(dat$time, fitted(gam_fit), lwd = 3) # Predictions

# points(dat$time, dat$y_obs)

# Mid left plot
gam_fit <- gam(y_obs ~ s(time, bs = "tp", k = 5), data = dat)
predmat <- predict(gam_fit, type = "lpmatrix")
predmat <- as.data.frame(t(t(predmat) * gam_fit$coefficients))

plot(dat$time, dat$y_obs,
  type = "n",
  main = "c) Five basis functions",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE
) # Data points

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

# Basis function values
for (i in seq_len(ncol(predmat))) {
  lines(
    x = dat$time, y = predmat[, i],
    col = rgb(red = 1, green = 0, blue = 0, alpha = 1)
  )
}
lines(dat$time, dat$y, lty = 2, lwd = 3) # State functions
lines(dat$time, fitted(gam_fit), lwd = 3) # Predictions

points(dat$time, dat$y_obs)

# Mid right plot
gam_fit <- gam(y_obs ~ s(time, bs = "tp", k = 100), data = dat)
predmat <- predict(gam_fit, type = "lpmatrix")
predmat <- as.data.frame(t(t(predmat) * gam_fit$coefficients))

plot(dat$time, dat$y_obs,
  type = "n",
  main = "d) One hundred basis functions",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE
) # Data points

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

# Basis function values
for (i in seq_len(ncol(predmat))) {
  lines(
    x = dat$time, y = predmat[, i],
    col = rgb(red = 1, green = 0, blue = 0, alpha = 1)
  )
}
lines(dat$time, dat$y, lty = 2, lwd = 3) # State functions
lines(dat$time, fitted(gam_fit), lwd = 3) # Predictions

points(dat$time, dat$y_obs)

# Bottom left plot
gam_fit <- gam(y_obs ~ s(time, bs = "tp", k = 20, sp = 0.0001), data = dat)
predmat <- predict(gam_fit, type = "lpmatrix")
predmat <- as.data.frame(t(t(predmat) * gam_fit$coefficients))

plot(dat$time, dat$y_obs,
  main = "e) Smoothing penalty too small",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE
) # Data points

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

# Basis function values
# for (i in seq_len(ncol(predmat))) {
#   lines(
#     x = dat$time, y = predmat[, i],
#     col = rgb(red = 1, green = 0, blue = 0, alpha = 1)
#   )
# }
lines(dat$time, dat$y, lty = 2, lwd = 3) # State functions
lines(dat$time, fitted(gam_fit), lwd = 3) # Predictions

# Bottom right plot
gam_fit <- gam(y_obs ~ s(time, bs = "tp", k = 20, sp = 1), data = dat)
predmat <- predict(gam_fit, type = "lpmatrix")
predmat <- as.data.frame(t(t(predmat) * gam_fit$coefficients))

plot(dat$time, dat$y_obs,
  main = "f) Smoothing penalty too large",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE
) # Data points

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

# Basis function values
# for (i in seq_len(ncol(predmat))) {
#   lines(
#     x = dat$time, y = predmat[, i],
#     col = rgb(red = 1, green = 0, blue = 0, alpha = 1)
#   )
# }
lines(dat$time, dat$y, lty = 2, lwd = 3) # State functions
lines(dat$time, fitted(gam_fit), lwd = 3) # Predictions

dev.off()

## GP -----
data <- list(
  N_obs = length(dat$time), x_obs = dat$time, y_obs = dat$y_obs,
  N_predict = length(dat$time),
  x_predict = dat$time
)

png(
  file = paste0(fig_path, "gp_demonstration.png"),
  width = 1960, height = 1080
)

layout_mat <- matrix(
  c(
    2, 2, 3, 3,
    4, 4, 5, 5,
    6, 6, 7, 7,
    1, 1, 1, 1
  ), 4, 4,
  byrow = TRUE
)
layout(mat = layout_mat, heights = c(0.3, 0.3, 0.3, 0.1))

par(lwd = 2.5)

# Top right plot
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
# Legend
legend(
  x = "bottom",
  legend = c("Process", "Prediction", "Sampled functions"),
  col = c("black", "black", "red"), lty = c(2, 1, 1), cex = cex,
  lwd = rep(2.5, 3),
  xpd = TRUE, horiz = TRUE, seg.len = 1, bty = "n"
)

par(lwd = 2.5)

# Top left plot
mod <- cmdstan_model("./R/stan_files/gp_prior.stan")

gp_fit <- mod$sample(
  data = data,
  fixed_param = TRUE,
  seed = 5838298,
  chains = 4,
  parallel_chains = 4,
  refresh = 500
)

# Posterior predictive draws
gp_post_pred(gp_fit,
  f_name = "f", time = dat$time,
  obs = dat$y_obs, state = rep(-20, 201), alpha = 0.2,
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex,
  main = "a) Gaussian process prior",
  xlab = "", ylab = "", axes = FALSE
)

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

lines(dat$time, dat$y, lty = 2, lwd = 3)
lines(x = dat$time, y = gp_fit$summary("f")$mean, lwd = 3)

# Top center plot
mod <- cmdstan_model("./R/stan_files/gp.stan")

gp_fit <- mod$sample(
  data = data,
  seed = 5838298,
  chains = 4,
  parallel_chains = 4,
  refresh = 500
)

# Posterior draws
gp_post_pred(gp_fit,
  f_name = "f_post_predict", time = dat$time,
  obs = dat$y_obs, state = rep(-20, 201), alpha = 0.2,
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex,
  main = "b) Gaussian process posterior",
  xlab = "", ylab = "", axes = FALSE
)

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

lines(dat$time, dat$y, lty = 2, lwd = 3)
lines(x = dat$time, y = gp_fit$summary("f_post_predict")$mean, lwd = 3)

# Mid left plot
mod <- cmdstan_model("./R/stan_files/gp_rho.stan")
data <- c(data, rho = 0.1)

gp_fit <- mod$sample(
  data = data,
  seed = 5838298,
  chains = 4,
  parallel_chains = 4,
  refresh = 500
)

# Posterior predictive draws
gp_post_pred(gp_fit,
  f_name = "f_post_predict", time = dat$time,
  obs = dat$y_obs, state = rep(-20, 201), alpha = 0.2,
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex,
  main = "c) Lengthscale too short",
  xlab = "", ylab = "", axes = FALSE
)

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

lines(dat$time, dat$y, lty = 2, lwd = 3)
lines(x = dat$time, y = gp_fit$summary("f_post_predict")$mean, lwd = 3)

# Mid right plot
mod <- cmdstan_model("./R/stan_files/gp_rho.stan")
data$rho <- 1

gp_fit <- mod$sample(
  data = data,
  seed = 5838298,
  chains = 4,
  parallel_chains = 4,
  refresh = 500
)

# Posterior predictive draws
gp_post_pred(gp_fit,
  f_name = "f_post_predict", time = dat$time,
  obs = dat$y_obs, state = rep(-20, 201), alpha = 0.2,
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex,
  main = "d) Lengthscale too large",
  xlab = "", ylab = "", axes = FALSE
)

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

lines(dat$time, dat$y, lty = 2, lwd = 3)
lines(x = dat$time, y = gp_fit$summary("f_post_predict")$mean, lwd = 3)

# Bottom left plot
mod <- cmdstan_model("./R/stan_files/gp_alpha.stan")
data <- c(data, alpha = 0.1)

gp_fit <- mod$sample(
  data = data,
  seed = 5838298,
  chains = 4,
  parallel_chains = 4,
  refresh = 500
)

# Posterior predictive draws
gp_post_pred(gp_fit,
  f_name = "f_post_predict", time = dat$time,
  obs = dat$y_obs, state = rep(-20, 201), alpha = 0.2,
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex,
  main = "e) Marginal variance too small",
  xlab = "", ylab = "", axes = FALSE
)

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

lines(dat$time, dat$y, lty = 2, lwd = 3)
lines(x = dat$time, y = gp_fit$summary("f_post_predict")$mean, lwd = 3)

# Bottom right plot
mod <- cmdstan_model("./R/stan_files/gp_alpha.stan")
data$alpha <- 30

gp_fit <- mod$sample(
  data = data,
  seed = 5838298,
  chains = 4,
  parallel_chains = 4,
  refresh = 500
)

# Posterior predictive draws
gp_post_pred(gp_fit,
  f_name = "f_post_predict", time = dat$time,
  obs = dat$y_obs, state = rep(-20, 201), alpha = 0.2,
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex,
  main = "f) Marginal variance too large",
  xlab = "", ylab = "", axes = FALSE
)

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

lines(dat$time, dat$y, lty = 2, lwd = 3)
lines(x = dat$time, y = gp_fit$summary("f_post_predict")$mean, lwd = 3)

dev.off()

# GP kernel plot
png(
  file = paste0(fig_path, "gp_kernel.png"),
  width = 1960, height = 1080
)

par(lwd = 2.5, mfrow = c(3, 1))

# Top plot
sqe <- function(x, alpha = 1, rho = 1) {
  alpha^2 * exp(-abs(x)^2 / (2 * rho^2))
}

loc <- seq(0, 4, 0.01)

plot(loc, sqe(loc),
  type = "l",
  main = "a) Squared exponential kernel",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE,
  lwd = 3
) # Data points

axis(2,
  at = c(0, 0.5, 1),
  labels = c("0", "", expression(alpha^2)), cex.axis = cex
)
axis(1,
  at = c(0, 1, 2, 3, 4),
  labels = c("0", expression(rho^2), "", expression(paste("3", rho^2)), ""),
  cex.axis = cex, padj = 1
)
title(xlab = expression(paste(Delta, "t")), line = 4, cex.lab = cex)

# Mid plot
matern_12 <- function(x, alpha = 1, rho = 1) {
  alpha^2 * exp(-abs(x / rho))
}

loc <- seq(0, 4, 0.01)

plot(loc, matern_12(loc),
  type = "l",
  main = "b) Matern 1/2 kernel",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE,
  lwd = 3
) # Data points

axis(2,
  at = c(0, 0.5, 1),
  labels = c("0", "", expression(alpha^2)), cex.axis = cex
)
axis(1,
  at = c(0, 1, 2, 3, 4),
  labels = c("0", expression(rho), "", expression(paste("3", rho)), ""),
  cex.axis = cex, padj = 1
)
title(xlab = expression(paste(Delta, "t")), line = 4, cex.lab = cex)


# Bottom plot
linear <- function(x, y, alpha = 1, beta = 0.01) {
  alpha^2 + beta^2 * x * y
}

matern_12 <- function(x, y, alpha = 1, rho = 50) {
  alpha^2 * exp(-abs(x - y) / rho)
}

sqe <- function(x, y, alpha = 1, rho = 50) {
  alpha^2 * exp(-abs(x - y)^2 / (2 * rho^2))
}

mu_fun <- function(x) {
  rep(0, length(x))
}

cov_fun <- function(x, y, kernel) {
  outer(x, y, kernel)
}

draws_list <- lapply(list(linear, matern_12, sqe), function(kernel, x) {
  as.matrix(mvrnorm(n = 3, mu_fun(x), cov_fun(x, x, kernel)))
}, x = seq(0, 200, 1))

draws_list <- do.call(rbind, draws_list)

plot(0,
  type = "n",
  main = "c) Prior draws from different kernels",
  xlab = "", ylab = "",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex, axes = FALSE,
  xlim = c(0, 200), ylim = c(min(draws_list), max(draws_list))
) # Data points

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

colors <- rep(c("black", "blue", "red"), each = 3)

for (i in seq_len(nrow(draws_list))) {
  lines(seq(0, 200, 1), draws_list[i, ], col = colors[i])
}

dev.off()

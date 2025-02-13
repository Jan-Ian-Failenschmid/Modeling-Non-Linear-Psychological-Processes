#' ----------------------------------------------------------------------------#
#' Title: Create Exemplar Plots                                                #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 25-03-2024                                                    #
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


### Set-up ---------------------------------------------------------------------
library(data.table)
# Load functions
invisible(sapply(
  c(paste0("./R/helper/", dir(path = "./R/helper"))),
  source
))
set.seed(41)
fig_path <- c("./figures/")
cex <- 3.5

### Create Exemplar Plot without Process Nosie ---------------------------------
png(
  file = paste0(fig_path, "/exemplar_no_process_noise.png"),
  width = 1960, height = 1080
)

par(mfrow = c(2, 2), lwd = 2.5)

# LCS ---
exp_growth <- new("gen_model",
  time = 200,
  model_name = "exp_growth",
  model_type = "DE",
  pars = list(yr = 0.02, ya = 2),
  delta = (50 / 7) / 234,
  stepsize = (50 / 7) / 6,
  start = list(
    formula(y ~ -2)
  ),
  state_eq = list(
    formula(y ~ yr * ya - yr * y)
  ),
  meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, 1)))
)

make_exemplar_plot(exp_growth,
  main = "a) Exponential Growth Curve",
  xlab = "", ylab = "", cex = cex, axes = FALSE
)
axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)

# Logistic growth ---
log_growth <- new("gen_model",
  time = 200,
  model_name = "log_growth",
  model_type = "DE",
  pars = list(k = 4.3, r = 0.04),
  delta = (50 / 7) / 234,
  stepsize = (50 / 7) / 6,
  start = list(
    formula(y ~ 0.3)
  ),
  state_eq = list(
    formula(y ~ r * y * (1 - (y / k)))
  ),
  meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, 1))),
  boundary = function(x) ifelse(x >= 0, x, 0)
)

make_exemplar_plot(log_growth,
  main = "b) Logistic Growth Curve",
  xlab = "", ylab = "", cex = cex, axes = FALSE
)

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)

# Cusp ---
cusp_catastrophe <- new("gen_model",
  time = 200,
  model_name = "catastrophe",
  model_type = "DE",
  delta = (50 / 7) / 234,
  stepsize = (50 / 7) / 6,
  pars = list(a = -5),
  start = list(
    formula(y ~ 1.9),
    formula(v ~ 0),
    formula(b ~ -10)
  ),
  state_eq = list(
    formula(y ~ -(4 * y^3 + 2 * a * y + b)),
    formula(b ~ v),
    formula(v ~ -(2 * pi / 50)^2 * b)
  ),
  meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, 1)))
)

make_exemplar_plot(cusp_catastrophe,
  main = "c) Cusp Catastrophe",
  xlab = "", ylab = "", cex = cex, axes = FALSE
)

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)

title(xlab = "Time", line = 4, cex.lab = cex)

# Dampened Oscillator ---
damped_oscillator <- new("gen_model",
  time = 200,
  model_name = "damped_oscillator",
  model_type = "DE",
  delta = (50 / 7) / 234,
  stepsize = (50 / 7) / 6,
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

make_exemplar_plot(damped_oscillator,
  main = "d) Damped Oscillator",
  xlab = "", ylab = "", cex = cex, axes = FALSE
)

axis(2, at = seq(-10, 10, 2), cex.axis = cex)
axis(1,
  at = c(0, 50, 100, 150, 200),
  labels = c("", "", "", "", ""),
  cex.axis = cex, padj = 1
)
title(xlab = "Time", line = 4, cex.lab = cex)

dev.off()

### Create Exemplar Plot with Process Nosie ------------------------------------
cex <- 4
png(
  file = paste0(fig_path, "exemplar_process_noise.png"),
  width = 1960, height = 1080
)

par(mfcol = c(4, 2))

for (dyn_er in sqrt(c(.5, 2))) {
  set.seed(1234)

  # LCS ---
  latent_change <- new("gen_model",
    time = 200,
    model_name = "lcs",
    model_type = "DE",
    pars = list(yr = 0.02, ya = 2, dyn_er = dyn_er),
    dynamic_error = list(
      formula(y ~ dyn_er)
    ),
    delta = (50 / 7) / 234,
    stepsize = (50 / 7) / 3,
    start = list(
      formula(y ~ -2)
    ),
    state_eq = list(
      formula(y ~ yr * ya - yr * y)
    ),
    meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, 1)))
  )

  make_exemplar_plot(latent_change,
    main = bquote(paste(
      "Exponential Growth Curve ", sigma^2, " = ",
      .(round(dyn_er^2, 1))
    )),
    xlab = "", ylab = "", cex = cex, axes = FALSE,
    col = rep(c("black", "blue"), each = 3 * 14),
    ylim = c(-10, 12)
  )

  axis(2, at = seq(-10, 20, 2), cex.axis = cex)
  axis(1,
    at = c(0, 50, 100, 150, 200),
    labels = c("", "", "", "", ""),
    cex.axis = cex, padj = 1
  )

  # Logistic growth ---
  log_growth <- new("gen_model",
    time = 200,
    model_name = "log_growth",
    model_type = "DE",
    pars = list(k = 4.3, r = 0.04, dyn_er = dyn_er),
    dynamic_error = list(
      formula(y ~ dyn_er)
    ),
    delta = (50 / 7) / 234,
    stepsize = (50 / 7) / 3,
    start = list(
      formula(y ~ 0.3)
    ),
    state_eq = list(
      formula(y ~ r * y * (1 - (y / k)))
    ),
    meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, 3))),
    boundary = function(x) ifelse(x >= 0, x, 0)
  )

  make_exemplar_plot(log_growth,
    main = bquote(paste(
      "Logistic Growth Curve ", sigma^2, " = ",
      .(round(dyn_er^2, 1))
    )),
    xlab = "", ylab = "", cex = cex, axes = FALSE,
    col = rep(c("black", "blue"), each = 3 * 14),
    ylim = c(-5, 14)
  )

  axis(2, at = seq(-10, 20, 2), cex.axis = cex)
  axis(1,
    at = c(0, 50, 100, 150, 200),
    labels = c("", "", "", "", ""),
    cex.axis = cex, padj = 1
  )

  # Cusp ----
  cusp_catastrophe <- new("gen_model",
    time = 200,
    model_name = "catastrophe",
    model_type = "DE",
    delta = (50 / 7) / 234,
    stepsize = (50 / 7) / 9,
    pars = list(a = -5, dyn_er = dyn_er),
    dynamic_error = list(
      formula(y ~ dyn_er),
      formula(b ~ 0),
      formula(v ~ 0)
    ),
    start = list(
      formula(y ~ 1.9),
      formula(v ~ 0),
      formula(b ~ -10)
    ),
    state_eq = list(
      formula(y ~ -(4 * y^3 + 2 * a * y + b)),
      formula(b ~ v),
      formula(v ~ -(2 * pi / 50)^2 * b)
    ),
    meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, 1)))
  )

  make_exemplar_plot(cusp_catastrophe,
    main = bquote(paste(
      "Cusp Catastrophe ", sigma^2, " = ",
      .(round(dyn_er^2, 1))
    )),
    xlab = "", ylab = "", cex = cex, axes = FALSE,
    col = rep(c("black", "blue"), each = 9 * 14),
    ylim = c(-5, 5)
  )

  axis(2, at = seq(-10, 10, 2), cex.axis = cex)
  axis(1,
    at = c(0, 50, 100, 150, 200),
    labels = c("", "", "", "", ""),
    cex.axis = cex, padj = 1
  )

  # Dampened Oscillator ---
  damp_osc <- new("gen_model",
    time = 200,
    model_name = "damped_oscillator",
    model_type = "DE",
    delta = (50 / 7) / 234,
    stepsize = (50 / 7) / 9,
    pars = list(b = 0.01, omega = 0.1, dyn_er = dyn_er),
    start = list(
      formula(y ~ 2),
      formula(v ~ 0)
    ),
    state_eq = list(
      formula(y ~ v),
      formula(v ~ -2 * b * v - (omega^2) * y)
    ),
    dynamic_error = list(
      formula(y ~ dyn_er),
      formula(v ~ 0)
    ),
    meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, 1)))
  )

  make_exemplar_plot(damp_osc,
    main = bquote(paste(
      "Damped Oscillator ", sigma^2, " = ",
      .(round(dyn_er^2, 1))
    )),
    xlab = "", ylab = "", cex = cex, axes = FALSE,
    col = rep(c("black", "blue"), each = 9 * 14),
    ylim = c(-16, 17)
  )

  axis(2, at = seq(-20, 20, 2), cex.axis = cex)
  axis(1,
    at = c(0, 50, 100, 150, 200),
    labels = c("", "First Half", "", "Second Half", ""),
    cex.axis = cex, padj = 1
  )
  title(xlab = "Time", line = 4, cex.lab = cex)
}

dev.off()

cat("Figures successfully generated!\n")

### Create summary exemplar plot with process noise ----------------------------

load("R/data/simulation_data_27_09_2024_00_32.Rdata")
sim[, model_name := sapply(gen_model, function(x) {
  x@model_name
})]

png(
  file = paste0(fig_path, "/stoch_mean_variance.png"),
  width = 1960, height = 1080
)

par(
  mfrow = c(4, 3), cex.lab = 3, cex.main = 3, cex.axis = 3,
  mar = c(5, 5, 3, 3) + .1
)
lay_mat <- matrix(c(
  3, 2, 1,
  6, 5, 4,
  9, 8, 7,
  12, 11, 10
), 4, 3, byrow = T)
layout(lay_mat)

for (model_name_iter in c(
  "exp_growth", "log_growth", "damped_oscillator", "cusp_catastrophe"
)) {
  for (dyn_er_iter in sqrt(c(2, 1, 0.5))) {
    dat <- sim[
      time == 200 &
        dyn_er == dyn_er_iter &
        stepsize == (50 / 7) / 9 &
        model_name == model_name_iter,
      dat
    ]

    x <- dat[[1]]$time

    dat <- do.call(rbind, lapply(dat, function(x) x$y))

    sum_stat <- t(apply(dat, 2, function(x) {
      c(
        quantile(x, probs = c((1 - 0.997) / 2, (1 - 0.95) / 2, (1 - 0.68) / 2)),
        mean(x),
        quantile(x, probs = c(
          1 - ((1 - 0.68) / 2), 1 - ((1 - 0.95) / 2), 1 - ((1 - 0.997) / 2)
        ))
      )
    }))

    if (dyn_er_iter == sqrt(2)) {
      ylim <- 1.1 * c(min(sum_stat), max(sum_stat))
    }

    plot(1,
      type = "n", xlab = "", ylab = "", axes = FALSE,
      xlim = c(0, 200), ylim = ylim
    )
    axis(2, at = seq(-20, 20, 2))
    axis(1,
      at = c(0, 50, 100, 150, 200),
      labels = c("", "First Half", "", "Second Half", ""),
      padj = 1
    )

    if (dyn_er_iter == sqrt(0.5)) {
      if (model_name_iter == "exp_growth") {
        title(main = expression(paste(sigma^2, " = 0.5")))
        title(ylab = "Exponential Growth")
      } else if (model_name_iter == "log_growth") {
        title(ylab = "Logistic Growth")
      } else if (model_name_iter == "damped_oscillator") {
        title(ylab = "Damped Oscillator")
      } else if (model_name_iter == "cusp_catastrophe") {
        title(ylab = "Cusp Catastrophe")
      }
    } else if (dyn_er_iter == sqrt(1) && model_name_iter == "exp_growth") {
      title(main = expression(paste(sigma^2, " = 1")))
    } else if (dyn_er_iter == sqrt(2) && model_name_iter == "exp_growth") {
      title(main = expression(paste(sigma^2, " = 2")))
    }

    if (model_name_iter == "cusp_catastrophe") {
      title(xlab = "Time", line = 4)
    }

    polygon(
      x = c(x, rev(x)), y = c(sum_stat[, 1], rev(sum_stat[, 7])),
      col = rgb(1, 0, 0, 0.3)
    )
    polygon(
      x = c(x, rev(x)), y = c(sum_stat[, 2], rev(sum_stat[, 6])),
      col = rgb(1, 0, 0, 0.6)
    )
    polygon(
      x = c(x, rev(x)), y = c(sum_stat[, 3], rev(sum_stat[, 5])),
      col = rgb(1, 0, 0, 1)
    )
    lines(x, sum_stat[, 4], lwd = 2.5)
  }
}

dev.off()

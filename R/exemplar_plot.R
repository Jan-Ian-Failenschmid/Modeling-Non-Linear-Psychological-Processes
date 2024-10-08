#' ----------------------------------------------------------------------------#
#' Title: Create Exemplar Plots                                                #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 25-03-2024                                                    #
#' -----                                                                       #
#' Last Modified: 08-10-2024                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#


### Set-up ---------------------------------------------------------------------
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

par(mfrow = c(2, 2))

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
    col = rep(c("black", "blue"), each = 3 * 14)
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
    col = rep(c("black", "blue"), each = 3 * 14)
  )

  axis(2, at = seq(-10, 10, 2), cex.axis = cex)
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
    col = rep(c("black", "blue"), each = 9 * 14)
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
    col = rep(c("black", "blue"), each = 9 * 14)
  )

  axis(2, at = seq(-10, 10, 2), cex.axis = cex)
  axis(1,
    at = c(0, 50, 100, 150, 200),
    labels = c("", "First Half", "", "Second Half", ""),
    cex.axis = cex, padj = 1
  )
  title(xlab = "Time", line = 4, cex.lab = cex)
}

dev.off()

cat("Figures successfully generated!\n")

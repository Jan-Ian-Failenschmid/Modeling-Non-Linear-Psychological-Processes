#' ----------------------------------------------------------------------------#
#' Title: Real Data Demonstration                                              #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 23-05-2024                                                    #
#' -----                                                                       #
#' Last Modified: 09-08-2024                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#

### Set-up ---------------------------------------------------------------------
# Load required packages
library(mgcv) # GAM's
library(cmdstanr) # Stan interface
library(dynr) # Dynamic modelling with regime switching
library(nprobust) # Local polynomial estimator
library(ggdist) # Plotting
library(data.table) # Data table for storing the simulation grid
library(marginaleffects) # GAMM probing
library(ggplot2) # Plotting
library(papaja) # Export package list

# Load functions
invisible(sapply(
  c(paste0("./R/helper/", dir(path = "./R/helper"))),
  source
))

### Functions ------------------------------------------------------------------
eval_points <- 200

plot_pred <- function(data, pred, ind, var) {
  time <- data[[ind]][, time]
  eval_time <- seq(min(time), max(time), length.out = eval_points)
  obs <- unlist((data[[ind]][, var, with = FALSE]))
  plot(x = time, y = obs)
  lines(x = eval_time, y = pred[[ind]]$est)
  lines(x = eval_time, y = pred[[ind]]$ub, lty = 2)
  lines(x = eval_time, y = pred[[ind]]$lb, lty = 2)
}

plot_ml_pred <- function(id, pred, dat, var, time) {
  pred_dat <- do.call(rbind.data.frame, mapply(function(UUID, pred) {
    data.frame(
      ID = UUID, time = pred$time, est = pred$est,
      ub = pred$ub, lb = pred$lb
    )
  }, UUID = id, pred = pred, SIMPLIFY = FALSE))

  dat_dat <- do.call(rbind.data.frame, mapply(
    function(UUID, dat, var, time) {
      data.frame(
        ID = UUID, val = dat[[var]], time = dat[[time]]
      )
    },
    UUID = id, dat = dat, MoreArgs = list(var = var, time = time),
    SIMPLIFY = FALSE
  ))

  gg <- ggplot(pred_dat, aes(x = time, y = est, col = ID, fill = ID)) +
    geom_line() +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = .1) +
    geom_point(data = dat_dat, mapping = aes(y = val, x = time)) +
    theme_apa() +
    theme(legend.position = "none")
  print(gg)
}

### Load in Data ---------------------------------------------------------------
# Read in data from secure storage location
inpd <- c("/mnt/c/Users/failensc/OneDrive - Tilburg University/Documenten/data")
inpf <- c("/data_downloads_XOVWITVGQ6_2023-11-09Leuven_clinical_study.csv")

df_raw <- fread(paste0(inpd, inpf))

df_raw[
  ,
  time := as.POSIXct(paste(
    as.Date(Date_Local, format = "%d/%m/%Y"),
    Time_Local
  ), tz = "CET")
]


df_raw[, time0 := difftime(time, min(time), units = "hours"), by = UUID]
df_raw[, time0 := difftime(time, min(time[weekdays(time) == "Monday"]),
  units = "hours"
), by = UUID]

df_raw$DEP_ES[df_raw$UUID == unique(df_raw$UUID)[110]]

df_raw$UUID <- as.factor(df_raw$UUID)

df <- df_raw[UUID != unique(UUID)[110], ] # Part 110 excluded due to no var

grouped_df <- df[, .(data = list(.SD)),
  by = UUID, .SDcols = c("DEP_ES", "time", "time0")
]

nrow(grouped_df) # Number of participants
summary(df[, .(age = unique(AGE_BL)), by = UUID][, age])
summary(df[, .(gender = unique(as.factor(GENDER_BL))), by = UUID][, gender])

### Fit Locpols ----------------------------------------------------------------
grouped_df[, locpol_fit := lapply(data, function(data) {
  # Remove missing data
  data <- data[!is.na(DEP_ES), ]
  # Find best polynomial degree
  # Evaluate kernel function with best degree
  loc_fit <- lprobust_cust(
    x = as.numeric(data$time0),
    y = data$DEP_ES,
    eval = as.numeric(seq(min(data[, time0]),
      max(data[, time0]),
      length.out = eval_points
    )),
    p = c(1), kernel = "gau", bwselect = "imse-dpi",
    bwcheck = 0, diag_A = FALSE
  )
})]

grouped_df[, locpol_pred := lapply(locpol_fit, function(fit) {
  inf <- as.data.frame(fit$Estimate)
  list(
    time = inf$eval,
    est = inf$tau.bc,
    ub = as.vector(inf$tau.bc + (qnorm(0.975) * inf$se.rb)),
    lb = as.vector(inf$tau.bc - (qnorm(0.975) * inf$se.rb))
  )
})]

grouped_df[, locpol_bw := lapply(locpol_fit, function(fit) {
  fit$Estimate[1, 2]
})]

grouped_df[, locpol_resid := mapply(
  function(data, locpol_pred) {
    data$DEP_ES - locpol_pred$est
  },
  data = data, locpol_pred = locpol_pred, SIMPLIFY = FALSE
)]

grouped_df[, locpol_mse := lapply(
  locpol_resid,
  function(resid) mean(resid^2, na.rm = TRUE)
)]

plot_pred(
  data = grouped_df[, data],
  pred = grouped_df[, locpol_pred], ind = 50, var = "DEP_ES"
) # 13, 22, 30

i <- 100:117
plot_ml_pred(
  id = grouped_df$UUID[i], pred = grouped_df$locpol_pred[i],
  dat = grouped_df$data[i], var = "DEP_ES", time = "time0"
)

which.min(grouped_df$locpol_bw)
### GP
grouped_df[, gp_fit := lapply(data, function(data) {
  data <- data[!is.na(DEP_ES), ]

  data$time0 <- (data$time0 - mean(data$time0)) / sd(data$time0)

  stan_data <- list(
    N = length(data$time0), x = as.numeric(data$time0), y = data$DEP_ES,
    N_eval = eval_points, x_eval = as.numeric(seq(min(data[, time0]),
      max(data[, time0]),
      length.out = eval_points
    )), c_f1 = 1.5,
    M_f1 = 50
  )

  mod <- quiet(cmdstan_model("./R/stan_files/gpbf1b.stan"))

  gp_fit <- suppressWarnings(mod$sample(
    data = stan_data,
    seed = 5838298,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    show_messages = TRUE,
    show_exceptions = FALSE
  ))
})]

grouped_df[, gp_pred := mapply(function(fit, data) {
  posterior_draws <- quiet(
    as.data.frame(fit$draws(
      "f",
      format = "draws_df"
    ))
  )

  posterior_draws <- posterior_draws[
    ,
    seq(1, ncol(posterior_draws) - 3)
  ]
  list(
    time = as.numeric(seq(min(data[, time0]),
      max(data[, time0]),
      length.out = eval_points
    )),
    est = sapply(posterior_draws, mean),
    ub = sapply(posterior_draws, quantile, probs = 0.975),
    lb = sapply(posterior_draws, quantile, probs = 0.025)
  )
}, fit = gp_fit, data = data, SIMPLIFY = FALSE)]

grouped_df[, gp_alpha := lapply(gp_fit, function(fit) {
  mean(fit$draws("sigma_f1", format = "draws_df")$sigma_f1)
})]

grouped_df[, gp_rho := lapply(gp_fit, function(fit) {
  mean(fit$draws("lengthscale_f1", format = "draws_df")$lengthscale_f1)
})]

grouped_df[, gp_sigma := lapply(gp_fit, function(fit) {
  mean(fit$draws("sigma", format = "draws_df")$sigma)
})]

grouped_df[, gp_resid := mapply(
  function(data, gp_pred) {
    data$DEP_ES - gp_pred$est
  },
  data = data, gp_pred = gp_pred, SIMPLIFY = FALSE
)]

grouped_df[, gp_mse := lapply(
  gp_resid,
  function(resid) mean(resid^2, na.rm = TRUE)
)]

# Compare participant 4
plot_pred(
  data = grouped_df[, data],
  pred = grouped_df[, gp_pred], ind = 6, var = "DEP_ES"
) # 13, 22, 30

plot_ml_pred(
  id = grouped_df$UUID[1:10], pred = grouped_df$gp_pred[1:10],
  dat = grouped_df$data[1:10], var = "DEP_ES", time = "time0"
)

cor(grouped_df$locpol_bw, grouped_df$gp_rho)
plot(x = grouped_df$locpol_bw, y = grouped_df$gp_rho)

### Fit GAMs -------------------------------------------------------------------
grouped_df[, gam_fit := lapply(data, function(data) {
  gam(
    DEP_ES ~ s(as.numeric(time0),
      bs = "tp",
      k = sum(!is.na(data$DEP_ES)) - 1
    ),
    data = data, method = "ML"
  )
})]

grouped_df[, gam_sp := lapply(gam_fit, function(fit) fit$sp)]

grouped_df[, gam_pred := mapply(function(fit, data) {
  eval_time <- data.frame(
    time0 = seq(min(data[, time0]), max(data[, time0]),
      length.out = eval_points
    )
  )
  inference <- predict(fit, newdata = eval_time, se.fit = TRUE)
  list(
    time = eval_time$time0,
    est = as.vector(inference$fit),
    ub = as.vector(inference$fit + (qnorm(0.975) * inference$se.fit)),
    lb = as.vector(inference$fit - (qnorm(0.975) * inference$se.fit))
  )
}, fit = gam_fit, data = data, SIMPLIFY = FALSE)]

grouped_df[, gam_resid := mapply(
  function(data, gam_pred) {
    data$DEP_ES - gam_pred$est
  },
  data = data, gam_pred = gam_pred, SIMPLIFY = FALSE
)]

grouped_df[, gam_mse := lapply(
  gam_resid,
  function(resid) mean(resid^2, na.rm = TRUE)
)]

plot_pred(
  data = grouped_df[, data],
  pred = grouped_df[, gam_pred], ind = 8, var = "DEP_ES"
) # 7 & 13

summary(grouped_df$gam_fit[[83]])

i <- c(52:63)
plot_ml_pred(
  id = grouped_df$UUID[i], pred = grouped_df$gam_pred[i],
  dat = grouped_df$data[i], var = "DEP_ES", time = "time0"
)

plot(x = grouped_df$locpol_bw[-108], y = grouped_df$gam_sp[-108])
summary(unlist(grouped_df$locpol_bw))
hist(unlist(grouped_df$gam_sp)[-108])
cor(unlist(grouped_df$locpol_bw), unlist(grouped_df$gam_sp))
order(unlist(grouped_df$gam_sp), decreasing = TRUE)
unlist(grouped_df$gam_sp)[108]
unlist(grouped_df$gam_sp)[83]
### GAM ------------------------------------------------------------------------
gamm <- gam(DEP_ES ~ s(as.numeric(time0)) +
  s(as.numeric(time0), UUID, bs = "fs"), data = df)
summary(gamm)
plot_predictions(gamm,
  condition = list(time0 = unique, UUID = unique(df$UUID)[1])
)


### Parametric modelling -------------------------------------------------------
set.seed(42)
alloc <- rep(
  sample(c("train", "valid", "test"), 117,
    prob = c(0.7, 0.15, 0.15), replace = TRUE
  ),
  each = 70
)

split_df <- split(df[, c("UUID", "DEP_ES", "time0")], alloc)
build <- split_df$train # [UUID == unique(UUID)[1], ]

# Random walk model -----
dynr_dat <- dynr.data(build, id = "UUID", time = "time0", observed = "DEP_ES")

meas <- prep.measurement(
  values.load = matrix(1, 1),
  params.load = matrix("fixed", 1),
  state.names = c("DEP"),
  obs.names = c("DEP_ES")
)

ecov <- prep.noise(
  values.latent = matrix(1, 1),
  params.latent = matrix("dyn_er", 1),
  values.observed = matrix(9, 1),
  params.observed = matrix("meas_er", 1)
)

dynamics <- prep.matrixDynamics(
  values.dyn = matrix(0, 1),
  params.dyn = matrix("fixed", 1),
  isContinuousTime = TRUE
)

initial <- prep.initial(
  values.inistate = na.omit(build$DEP_ES)[[1]],
  params.inistate = c("DEP_nod"),
  values.inicov = matrix(9, 1),
  params.inicov = matrix("DEP_nod_var", 1)
)

model <- dynr.model(
  dynamics = dynamics, measurement = meas, noise = ecov,
  initial = initial, data = dynr_dat
)

model$lb[
  c("dyn_er", "meas_er", "DEP_nod")
] <- c(
  1e-6, 1e-6, 0
)

model$ub[
  c("dyn_er", "meas_er", "DEP_nod")
] <- c(
  1e+3, 1e+3, 1e+2
)

dynr_fit <- dynr.cook(model)
summary(dynr_fit)
plot(dynr_fit, model, names.observed = "DEP_ES")
coef(model) <- coef(dynr_fit)
prediction <- predict(model)$estimate
build$DEP_hat <- prediction[1, ]

ggplot(data = build, aes(y = DEP_ES, x = time0, col = UUID)) +
  geom_point() +
  geom_line(aes(y = DEP_hat)) +
  theme_apa() +
  theme(legend.position = "none")

# AR(1) -----
dynr_dat <- dynr.data(build, id = "UUID", time = "time0", observed = "DEP_ES")

meas <- prep.measurement(
  values.load = matrix(c(1), 1, 1),
  params.load = matrix("fixed", 1, 1),
  state.names = c("DEP"),
  obs.names = c("DEP_ES")
)

ecov <- prep.noise(
  values.latent = diag(c(1), 1),
  params.latent = diag(c("dyn_er"), 1),
  values.observed = matrix(9, 1),
  params.observed = matrix("meas_er", 1)
)

dynamics <- prep.formulaDynamics(
  formula = list(
    DEP ~ rho * (mu - DEP)
  ),
  startval = c(mu = 50, rho = 0.5),
  isContinuousTime = TRUE
  # isContinuousTime = TRUE,
  # theta.formula = list(mu ~ rand_par),
  # random.names = c("rand_par"),
  # random.params.inicov = matrix("mu_var", 1),
  # random.values.inicov = matrix(9, 1)
)

initial <- prep.initial(
  values.inistate = c(50),
  params.inistate = c("DEP_nod"),
  values.inicov = diag(c(10), 1),
  params.inicov = diag(c("DEP_nod_var"), 1)
)

model <- dynr.model(
  dynamics = dynamics, measurement = meas, noise = ecov,
  initial = initial, data = dynr_dat
)

model$lb[
  c("mu", "rho", "dyn_er", "meas_er", "DEP_nod", "DEP_nod_var")
] <- c(
  0, 0, 1e-6, 1e-6, 0, 1e-6
)

model$ub[
  c("mu", "rho", "dyn_er", "meas_er", "DEP_nod", "DEP_nod_var")
] <- c(
  1e+2, 1, 1e+3, 1e+3, 1e+2, 1e+3
)

dynr_fit <- dynr.cook(model)
summary(dynr_fit)
coef(model) <- coef(dynr_fit)
prediction <- predict(model)$estimate
build$DEP_hat <- prediction[1, ]
build$mu_hat <- prediction[2, ]

ggplot(data = build, aes(y = DEP_ES, x = time0, col = UUID)) +
  geom_point() +
  geom_line(aes(y = mu_hat)) +
  theme_apa() +
  theme(legend.position = "none")

# Damped Oscillator -----
meas <- prep.measurement(
  values.load = matrix(c(1, 0), 1, 2),
  params.load = matrix(c("fixed", "fixed"), 1, 2),
  state.names = c("DEP", "val"),
  obs.names = c("DEP_ES")
)

ecov <- prep.noise(
  values.latent = diag(c(1, 0), 2),
  params.latent = diag(c("dyn_er", "fixed"), 2),
  values.observed = diag(50, 1),
  params.observed = diag("meas_er", 1)
)

dynamics <- prep.formulaDynamics(
  formula = list(
    DEP ~ val,
    val ~ k * (DEP - mu) + c * val
  ),
  startval = c(mu = 50, k = 0, c = 0),
  isContinuousTime = TRUE
)

initial <- prep.initial(
  values.inistate = c(na.omit(build$DEP_ES)[[1]], 0),
  params.inistate = c("DEP_nod", "val_nod"),
  values.inicov = diag(c(1, 1), 2),
  params.inicov = diag(c("DEP_nod_var", "val_nod_var"), 2)
)

model <- dynr.model(
  dynamics = dynamics, measurement = meas, noise = ecov,
  initial = initial, data = dynr_dat
)

model$lb[
  c("k", "c", "mu", "dyn_er", "meas_er", "DEP_nod", "val_nod")
] <- c(
  -1e2, -1e2, 0, 1e-6, 1e-6, 0, -1e3
)

model$ub[
  c("k", "c", "mu", "dyn_er", "meas_er", "DEP_nod", "val_nod")
] <- c(
  0, 0, 1e+2, 1e+3, 1e+3, 1e+2, 1e3
)

dynr_fit <- dynr.cook(model)
summary(dynr_fit)

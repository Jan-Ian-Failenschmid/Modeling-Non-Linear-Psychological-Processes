#' ----------------------------------------------------------------------------#
#' Title: Real Data Demonstration                                              #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 23-05-2024                                                    #
#' -----                                                                       #
#' Last Modified: 21-08-2024                                                   #
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
# install.packages("mgcv") # GAM's
# install.packages("cmdstanr") # Stan interface
# install.packages("invgamma")
# install.packages("dynr") # Dynamic modelling with regime switching
# install.packages("nprobust") # Local polynomial estimator
# install.packages("ggdist") # Plotting
# install.packages("data.table") # Data table for storing the simulation grid
# install.packages("marginaleffects") # GAMM probing
# install.packages("ggplot2") # Plotting
# install.packages("patchwork")
# install.packages("papaja") # Export package list

library(mgcv) # GAM's
library(cmdstanr) # Stan interface
library(dynr) # Dynamic modelling with regime switching
library(nprobust) # Local polynomial estimator
library(ggdist) # Plotting
library(data.table) # Data table for storing the simulation grid
library(marginaleffects) # GAMM probing
library(ggplot2) # Plotting
library(patchwork)
library(papaja) # Export package list

# Load functions
invisible(sapply(
  c(paste0("./R/helper/", dir(path = "./R/helper"))),
  source
))

### Functions ------------------------------------------------------------------
plot_pred <- function(data, pred, ind, var) {
  data <- data[[ind]][!is.na(DEP_ES), ]
  time <- data[, time]
  if ("eval_points" %in% ls()) {
    eval_time <- seq(min(time), max(time), length.out = eval_points)
  } else {
    eval_time <- time
  }
  obs <- unlist((data[, var, with = FALSE]))
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
    theme(text = element_text(size = 20))
  return(gg)
}

plot_points <- 300

### Load in Data ---------------------------------------------------------------
# Read in data from secure storage location
inpd <- c("/mnt/c/Users/failensc/OneDrive - Tilburg University/Documenten/data")
inpd <- c("/mnt/c/Users/janfa/OneDrive - Tilburg University/Documenten/data")
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

min_monday <- min(df_raw$time[weekdays(df_raw$time) == "Monday"])
df_raw[, time1 := difftime(time, min_monday,
  units = "hours"
), by = UUID]
df_raw[, time1 := time1 -
  ceiling(min(time1) / (24 * 7)) * (24 * 7), by = UUID]

df_raw$UUID <- as.factor(df_raw$UUID)
levels(df_raw$UUID) <- 1:118
df <- df_raw[UUID != unique(UUID)[110], ] # Part 110 excluded due to no var

grouped_df <- df[, .(data = list(.SD)),
  by = UUID, .SDcols = c("DEP_ES", "time", "time0")
]

# Descriptives
nrow(grouped_df) # Number of participants
summary(df[, .(age = unique(AGE_BL)), by = UUID][, age])
summary(df[, .(gender = unique(as.factor(GENDER_BL))), by = UUID][, gender])


#### Idiographic modelling
### Fit Locpols ----------------------------------------------------------------
grouped_df[, locpol_fit := lapply(data, function(data) {
  # Remove missing data
  data <- data[!is.na(DEP_ES), ]
  # Find best polynomial degree
  # Evaluate kernel function with best degree
  loc_fit <- lprobust_cust(
    x = as.numeric(data$time0),
    y = data$DEP_ES,
    eval = as.numeric(data$time0),
    p = 3, kernel = "gau", bwselect = "imse-dpi",
    bwcheck = 0, diag_A = TRUE
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
    na.omit(data$DEP_ES) - locpol_pred$est
  },
  data = data, locpol_pred = locpol_pred, SIMPLIFY = FALSE
)]

grouped_df[, locpol_mse := lapply(
  locpol_resid,
  function(resid) mean(resid^2, na.rm = TRUE)
)]

grouped_df[, locpol_gcv := mapply(
  function(fit, data) {
    get_cv(fit$Estimate, na.omit(data$DEP_ES))
  },
  fit = locpol_fit, data = data, SIMPLIFY = FALSE
)]

grouped_df[, locpol_fit_plot := lapply(data, function(data) {
  # Remove missing data
  data <- data[!is.na(DEP_ES), ]
  # Find best polynomial degree
  # Evaluate kernel function with best degree
  loc_fit <- lprobust_cust(
    x = as.numeric(data$time0),
    y = data$DEP_ES,
    eval = seq(min(data$time0), max(data$time0), length.out = plot_points),
    p = 3, kernel = "gau", bwselect = "imse-dpi",
    bwcheck = 0, diag_A = FALSE
  )
})]

grouped_df[, locpol_pred_plot := mapply(function(fit, data) {
  data <- data[!is.na(DEP_ES), ]
  inf <- as.data.frame(fit$Estimate)
  list(
    time = seq(min(data$time0), max(data$time0), length.out = plot_points),
    est = inf$tau.bc,
    ub = as.vector(inf$tau.bc + (qnorm(0.975) * inf$se.rb)),
    lb = as.vector(inf$tau.bc - (qnorm(0.975) * inf$se.rb))
  )
}, fit = locpol_fit_plot, data = grouped_df$data, SIMPLIFY = FALSE)]

### GP
grouped_df[, gp_fit := lapply(data, function(data) {
  data <- data[!is.na(DEP_ES), ]

  # t_diff <- as.numeric(abs(outer(scale(data$time0), scale(data$time0),
  #   FUN = "-"
  # )[
  #   lower.tri(outer(data$time0, data$time0, FUN = "-"))
  # ]))

  stan_data <- list(
    N_obs = length(data$time0), x_obs = as.numeric(data$time0),
    y_obs = data$DEP_ES # , t_diff_min = 2, t_diff_max = max(t_diff)
  )

  mod <- quiet(cmdstan_model("./R/stan_files/gp.stan"))

  gp_fit <- mod$sample(
    data = stan_data,
    seed = 5838298,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    show_messages = TRUE,
    show_exceptions = TRUE
  )

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

  gp_pred <- list(
    time = data$time0,
    est = sapply(posterior_draws, mean),
    ub = sapply(posterior_draws, quantile, probs = 0.975),
    lb = sapply(posterior_draws, quantile, probs = 0.025)
  )

  gp_alpha <- mean(gp_fit$draws("alpha", format = "draws_df")$alpha)
  gp_rho <- mean(gp_fit$draws("rho", format = "draws_df")$rho)
  gp_sigma <- mean(gp_fit$draws("sigma", format = "draws_df")$sigma)
  gp_gcv <- mean(gp_fit$draws("gcv_val", format = "draws_df")$gcv_val)

  return(list(gp_pred, gp_alpha, gp_rho, gp_sigma, gp_gcv))
})]

grouped_df[, gp_pred := lapply(gp_fit, "[[", 1)]
grouped_df[, gp_alpha := lapply(gp_fit, "[[", 2)]
grouped_df[, gp_rho := lapply(gp_fit, "[[", 3)]
grouped_df[, gp_sigma := lapply(gp_fit, "[[", 4)]
grouped_df[, gp_gcv := lapply(gp_fit, "[[", 5)]

grouped_df[, gp_resid := mapply(
  function(data, gp_pred) {
    na.omit(data$DEP_ES) - gp_pred$est
  },
  data = data, gp_pred = gp_pred, SIMPLIFY = FALSE
)]

grouped_df[, gp_mse := lapply(
  gp_resid,
  function(resid) mean(resid^2, na.rm = TRUE)
)]

grouped_df[, gp_fit_plot := lapply(data, function(data) {
  data <- data[!is.na(DEP_ES), ]

  # t_diff <- as.numeric(abs(outer(scale(data$time0), scale(data$time0),
  #   FUN = "-"
  # )[
  #   lower.tri(outer(data$time0, data$time0, FUN = "-"))
  # ]))

  stan_data <- list(
    N_obs = length(data$time0), x_obs = as.numeric(data$time0),
    y_obs = data$DEP_ES, N_pred = plot_points,
    x_pred = seq(min(data$time0), max(data$time0), length.out = plot_points)
    # , t_diff_min = 2, t_diff_max = max(t_diff)
  )

  mod <- quiet(cmdstan_model("./R/stan_files/gp_plot.stan"))

  gp_fit <- mod$sample(
    data = stan_data,
    seed = 5838298,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    show_messages = TRUE,
    show_exceptions = TRUE
  )

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

  gp_pred <- list(
    time = seq(min(data$time0), max(data$time0), length.out = plot_points),
    est = sapply(posterior_draws, mean),
    ub = sapply(posterior_draws, quantile, probs = 0.975),
    lb = sapply(posterior_draws, quantile, probs = 0.025)
  )

  return(list(gp_pred))
})]

grouped_df[, gp_pred_plot := lapply(gp_fit_plot, "[[", 1)]

### Fit GAMs -------------------------------------------------------------------
grouped_df[, gam_fit := lapply(data, function(data) {
  data <- data[!is.na(DEP_ES), ]
  gam(
    DEP_ES ~ s(as.numeric(time0),
      bs = "tp",
      k = length(data$DEP_ES) - 1
    ),
    data = data, method = "ML"
  )
})]

grouped_df[, gam_pred := mapply(function(fit, data) {
  inference <- predict(fit, se.fit = TRUE)
  data <- data[!is.na(DEP_ES), ]
  list(
    time = data$time0,
    est = as.vector(inference$fit),
    ub = as.vector(inference$fit + (qnorm(0.975) * inference$se.fit)),
    lb = as.vector(inference$fit - (qnorm(0.975) * inference$se.fit))
  )
}, fit = gam_fit, data = data, SIMPLIFY = FALSE)]

grouped_df[, gam_sp := lapply(gam_fit, function(fit) {
  fit$sp
})]

grouped_df[, gam_resid := mapply(
  function(data, gam_pred) {
    na.omit(data$DEP_ES) - gam_pred$est
  },
  data = data, gam_pred = gam_pred, SIMPLIFY = FALSE
)]

grouped_df[, gam_mse := lapply(
  gam_resid,
  function(resid) mean(resid^2, na.rm = TRUE)
)]

grouped_df[, gam_gcv := mapply(
  function(fit, resid) {
    n <- length(resid)
    n * sum(resid^2) / (n - sum(influence(fit)))^2
  },
  fit = gam_fit, resid = gam_resid
)]

grouped_df[, gam_sigma := lapply(
  gam_fit,
  function(fit) {
    sqrt(sum(fit$residuals^2) / fit$df.residual)
  }
)]

grouped_df[, gam_pred_plot := mapply(function(fit, data) {
  data <- data[!is.na(DEP_ES), ]
  inference <- predict(fit,
    newdata = data.frame(
      time0 = seq(min(data$time0), max(data$time0), length.out = plot_points)
    ),
    se.fit = TRUE
  )
  list(
    time = seq(min(data$time0), max(data$time0), length.out = plot_points),
    est = as.vector(inference$fit),
    ub = as.vector(inference$fit + (qnorm(0.975) * inference$se.fit)),
    lb = as.vector(inference$fit - (qnorm(0.975) * inference$se.fit))
  )
}, fit = gam_fit, data = data, SIMPLIFY = FALSE)]

### Results --------------------------------------------------------------------
res_df <- data.table(
  method = as.factor(c(
    rep("locpol", nrow(grouped_df)),
    rep("gp", nrow(grouped_df)),
    rep("gam", nrow(grouped_df))
  )),
  wiggliness = c(
    unlist(grouped_df$locpol_bw),
    unlist(grouped_df$gp_rho),
    unlist(grouped_df$gam_sp)
  ),
  mse = c(
    unlist(grouped_df$locpol_mse),
    unlist(grouped_df$gp_mse),
    unlist(grouped_df$gam_mse)
  ),
  gcv = c(
    unlist(grouped_df$locpol_gcv),
    unlist(grouped_df$gp_gcv),
    unlist(grouped_df$gam_gcv)
  )
)

res_df[, lapply(.SD, function(x) round(mean(x), 2)), by = method]
res_df[, lapply(.SD, function(x) round(sd(x), 2)), by = method]
res_df[, lapply(.SD, function(x) round(median(x), 2)), by = method]
res_df[, lapply(.SD, function(x) round(IQR(x), 2)), by = method]

data.table(
  bandwidth = unlist(grouped_df$locpol_bw),
  lengthscale = unlist(grouped_df$gp_rho),
  smoothing_par = unlist(grouped_df$gam_sp)
) |> cor()

np <- 10
locpol_high <- plot_ml_pred(
  id = grouped_df$UUID[order(unlist(grouped_df$locpol_bw))[1:np]],
  pred = grouped_df$locpol_pred_plot[order(unlist(grouped_df$locpol_bw))[1:np]],
  dat = grouped_df$data[order(unlist(grouped_df$locpol_bw))[1:np]],
  var = "DEP_ES", time = "time0"
)
locpol_high <- locpol_high +
  ggtitle("LPRs with the lowest bandwidth") +
  xlab("Momentary depression") +
  ylab("Differenced time")

locpol_low <- plot_ml_pred(
  id = grouped_df$UUID[order(unlist(grouped_df$locpol_bw),
    decreasing = TRUE
  )[1:np]],
  pred = grouped_df$locpol_pred_plot[order(unlist(grouped_df$locpol_bw),
    decreasing = TRUE
  )[1:np]],
  dat = grouped_df$data[order(unlist(grouped_df$locpol_bw),
    decreasing = TRUE
  )[1:np]],
  var = "DEP_ES", time = "time0"
)
locpol_low <- locpol_low +
  ggtitle("LPRs with the highest bandwidth") +
  xlab("Momentary depression") +
  ylab("Differenced time")


gp_high <- plot_ml_pred(
  id = grouped_df$UUID[order(unlist(grouped_df$gp_rho))[1:np]],
  pred = grouped_df$gp_pred_plot[order(unlist(grouped_df$gp_rho))[1:np]],
  dat = grouped_df$data[order(unlist(grouped_df$gp_rho))[1:np]],
  var = "DEP_ES", time = "time0"
)
gp_high <- gp_high +
  ggtitle("GPs with the lowest lengthscale") +
  xlab("Momentary depression") +
  ylab("Differenced time")


gp_low <- plot_ml_pred(
  id = grouped_df$UUID[order(unlist(grouped_df$gp_rho),
    decreasing = TRUE
  )[1:np]],
  pred = grouped_df$gp_pred_plot[order(unlist(grouped_df$gp_rho),
    decreasing = TRUE
  )[1:np]],
  dat = grouped_df$data[order(unlist(grouped_df$gp_rho),
    decreasing = TRUE
  )[1:np]],
  var = "DEP_ES", time = "time0"
)
gp_low <- gp_low +
  ggtitle("GPs with the highest lengthscale") +
  xlab("Momentary depression") +
  ylab("Differenced time")

gam_high <- plot_ml_pred(
  id = grouped_df$UUID[order(unlist(grouped_df$gam_sp))[1:np]],
  pred = grouped_df$gam_pred_plot[order(unlist(grouped_df$gam_sp))[1:np]],
  dat = grouped_df$data[order(unlist(grouped_df$gam_sp))[1:np]],
  var = "DEP_ES", time = "time0"
)
gam_high <- gam_high +
  ggtitle("GAMs with the lowest smoothing parameter") +
  xlab("Momentary depression") +
  ylab("Differenced time")

gam_low <- plot_ml_pred(
  id = grouped_df$UUID[order(unlist(grouped_df$gam_sp),
    decreasing = TRUE
  )[1:np]],
  pred = grouped_df$gam_pred_plot[order(unlist(grouped_df$gam_sp),
    decreasing = TRUE
  )[1:np]],
  dat = grouped_df$data[order(unlist(grouped_df$gam_sp),
    decreasing = TRUE
  )[1:np]],
  var = "DEP_ES", time = "time0"
)
gam_low <- gam_low +
  ggtitle("GAMs with the highest smoothing parameter") +
  xlab("Momentary depression") +
  ylab("Differenced time")

complete_plot <- (locpol_low | gp_low | gam_low) /
  (locpol_high | gp_high | gam_high)

ggsave("figures/demonstration_smooths.png", complete_plot,
  width = 1920,
  height = 1080, units = "px", dpi = "screen"
)


#### Multilevel modelling
### GAMM------------------------------------------------------------------------
gamm_ri <- gam(DEP_ES ~ s(UUID, bs = "re"), data = df)
gamm_ncs <- gam(DEP_ES ~ s(as.numeric(time1), UUID, bs = "fs"), data = df)
gamm <- gam(DEP_ES ~ s(as.numeric(time1)) +
  s(as.numeric(time1), UUID, bs = "fs"), data = df)
gamm_free <- gam(DEP_ES ~ s(as.numeric(time1), by = UUID), data = df)

summary(gamm_ri)
summary(gamm_ncs)
summary(gamm)
summary(gamm_free)

AIC(gamm_ri, gamm_ncs, gamm, gamm_free)
BIC(gamm_ri, gamm_ncs, gamm, gamm_free)

plot(gamm, select = 2)
plot_predictions(gamm,
  condition = list(time0 = unique, UUID = unique(df$UUID)[1])
)

#### Parametric modelling ------------------------------------------------------
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

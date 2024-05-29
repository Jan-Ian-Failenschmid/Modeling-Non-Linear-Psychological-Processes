#' ----------------------------------------------------------------------------#
#' Title: Real Data Demonstration                                              #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 23-05-2024                                                    #
#' -----                                                                       #
#' Last Modified: 24-05-2024                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#

### Set-up ---------------------------------------------------------------------
if (!require(pacman)) install.packages("pacman")
pacman::p_load(
  mgcv, # GAM's
  cmdstanr, # Stan interface
  dynr, # Dynamic
  nprobust, # Local polynomial estimator
  data.table, # Data table for storing the simulation grid
  ggdist # Plotting
)

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

### Load in Data ---------------------------------------------------------------
# Read in data from secure storage location
inpd <- c("/mnt/c/Users/failensc/OneDrive - Tilburg University/Documenten/data")
inpf <- c("/data_downloads_XOVWITVGQ6_2023-11-09Leuven_clinical_study.csv")
df <- fread(paste0(inpd, inpf))

df[
  ,
  time := as.POSIXct(paste(
    as.Date(Date_Local, format = "%d/%m/%Y"),
    Time_Local
  ), tz = "CET")
]

grouped_df <- df[, .(data = list(.SD)),
  by = UUID, .SDcols = c("DEP_ES", "time")
]
grouped_df <- grouped_df[-110, ]

### Fit GAMs -------------------------------------------------------------------
grouped_df[, gam_fit := lapply(data, function(data) {
  gam(DEP_ES ~ s(as.numeric(time),
    bs = "tp", k = 34
  ), data = data)
})]

grouped_df[, gam_pred := mapply(function(fit, data) {
  eval_time <- data.frame(
    time = seq(min(data[, time]), max(data[, time]), length.out = eval_points)
  )
  inference <- predict(fit, newdata = eval_time, se.fit = TRUE)
  list(
    est = as.vector(inference$fit),
    ub = as.vector(inference$fit + (qnorm(0.975) * inference$se.fit)),
    lb = as.vector(inference$fit - (qnorm(0.975) * inference$se.fit))
  )
}, fit = gam_fit, data = data, SIMPLIFY = FALSE)]

plot_pred(
  data = grouped_df[, data],
  pred = grouped_df[, gam_pred], ind = 2, var = "DEP_ES"
) # 13, 22, 30

### Fit Locpols ----------------------------------------------------------------
grouped_df[, locpol_fit := lapply(data, function(data) {
  # Remove missing data
  data <- data[!is.na(DEP_ES), ]
  # Find best polynomial degree
  loc_fit_list <- lapply(c(1, 3, 5), function(p, method, data) {
    loc_fit <- lprobust_cust(
      x = as.numeric(data$time), y = data$DEP_ES, eval = as.numeric(data$time),
      p = p, kernel = "gau", bwselect = "imse-dpi",
      bwcheck = 0, diag_A = TRUE
    )
    gcv <- get_cv(loc_fit$Estimate, data$y_obs)
    return(list(loc_fit, gcv))
  }, data = data)
  ind <- which.min(sapply(loc_fit_list, "[[", 2))
  # Evaluate kernel function with best degree
  loc_fit <- lprobust_cust(
    x = as.numeric(data$time),
    y = data$DEP_ES,
    eval = as.numeric(seq(min(data[, time]),
      max(data[, time]),
      length.out = eval_points
    )),
    p = c(1, 3, 5)[ind], kernel = "gau", bwselect = "imse-dpi",
    bwcheck = 0, diag_A = FALSE
  )
})]

grouped_df[, locpol_pred := lapply(locpol_fit, function(fit) {
  inf <- as.data.frame(fit$Estimate)
  list(
    est = inf$tau.bc,
    ub = as.vector(inf$tau.bc + (qnorm(0.975) * inf$se.rb)),
    lb = as.vector(inf$tau.bc - (qnorm(0.975) * inf$se.rb))
  )
})]

plot_pred(
  data = grouped_df[, data],
  pred = grouped_df[, locpol_pred], ind = 2, var = "DEP_ES"
) # 13, 22, 30

### GP
grouped_df[, gp_fit := lapply(data, function(data) {
  data <- data[!is.na(DEP_ES), ]

  data$time <- (data$time - mean(data$time)) / sd(data$time)

  stan_data <- list(
    N = length(data$time), x = as.numeric(data$time), y = data$DEP_ES,
    N_eval = eval_points, x_eval = as.numeric(seq(min(data[, time]),
      max(data[, time]),
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

grouped_df[, gp_pred := lapply(gp_fit, function(fit) {
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
    est = sapply(posterior_draws, mean),
    ub = sapply(posterior_draws, quantile, probs = 0.975),
    lb = sapply(posterior_draws, quantile, probs = 0.025)
  )
})]

plot_pred(
  data = grouped_df[, data],
  pred = grouped_df[, gam_pred], ind = 1, var = "DEP_ES"
) # 13, 22, 30

### Analyze one data set -------------------------------------------------------
data <- df[UUID == UUID[1], ]

gam_fit <- gam(DEP_ES ~ s(as.numeric(time), bs = "tp", k = 34), data = data)
summary(gam_fit)
plot(gam_fit, residuals = TRUE)
acf(resid(gam_fit), lag.max = 10, main = "ACF")

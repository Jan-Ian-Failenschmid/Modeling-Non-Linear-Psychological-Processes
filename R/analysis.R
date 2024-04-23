#' ----------------------------------------------------------------------------#
#' Title: Main Analyisis                                                       #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 12-04-2024                                                    #
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

### Dependencies ---------------------------------------------------------------
if (!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, ggdist)

# Load functions
invisible(sapply(
  c(paste0("./R/helper/", dir(path = "./R/helper"))),
  source
))

### Load data ------------------------------------------------------------------
# load("R/data/simulation_results_23_04_2024_12_26.Rdata")

### Analysis -------------------------------------------------------------------
## Descriptives
res_summary <- res %>%
  group_by(method, model, time, stepsize, dyn_er) %>%
  summarise(
    mean_mse = mean(mse, na.rm = TRUE),
    sd_mse = sd(mse, na.rm = TRUE),
    mean_gcv = mean(gcv, na.rm = TRUE),
    sd_gcv = sd(gcv, na.rm = TRUE),
    mean_ci_coverage = mean(ci_coverage, na.rm = TRUE),
    sd_ci_coverage = sd(ci_coverage, na.rm = TRUE)
  )

View(res_summary)

## Visulization
x11()
plot_results(res = res, "mse", time, stepsize)
plot_results(res = res, "gcv", time, stepsize)
plot_results(res = res, "ci_coverage", time, stepsize)

## Inference
predictor <- c("method", "model", "time", "stepsize", "dyn_er")

# Model selection - Univariate
# mse
fit_fun_mse <- function(rhs, data) {
  form <- as.formula(paste0("mse ~ ", rhs))
  fit <- lm(form, data = data)
  return(fit)
}

mse_sel <- model_search(fit_fun_mse, predictor, data = res)
plot(seq_len(nrow(mse_sel)), mse_sel$w_AIC, type = "l")
plot(seq_len(nrow(mse_sel)), mse_sel$w_BIC, type = "l")

# GCV
fit_fun_gcv <- function(rhs, data) {
  form <- as.formula(paste0("gcv ~ ", rhs))
  fit <- lm(form, data = data)
  return(fit)
}

gcv_sel <- model_search(fit_fun_gcv, predictor, data = res)
plot(seq_len(nrow(gcv_sel)), gcv_sel$w_AIC, type = "l")
plot(seq_len(nrow(gcv_sel)), gcv_sel$w_BIC, type = "l")

# Model fitting
# AIC
fit_aic_mse <- lm(
  as.formula(
    paste0("mse ~ ", mse_sel$predictors[which.max(m_aic)])
  ),
  data = res
)
anova(fit_aic_mse)

fit_aic_gcv <- lm(
  as.formula(
    paste0("gcv ~ ", mse_sel$predictors[which.max(m_aic)])
  ),
  data = res
)
anova(fit_aic_gcv)

# BIC
fit_bic <- manova(
  as.formula(
    paste0("cbind(mse, gcv) ~ ", mse_sel$predictors[which.max(m_bic)])
  ),
  data = res
)
summary(fit_bic)

fit_bic_mse <- lm(
  as.formula(
    paste0("mse ~ ", mse_sel$predictors[which.max(m_bic)])
  ),
  data = res
)
anova(fit_bic_mse)

fit_bic_gcv <- lm(
  as.formula(
    paste0("gcv ~ ", mse_sel$predictors[which.max(m_bic)])
  ),
  data = res
)
anova(fit_bic_gcv)

### GCV vs. OCV
log_growth <- new("gen_model",
  model_name = "log_growth",
  model_type = "DE",
  time = 200,
  stepsize = 1,
  pars = list(k = 4.3, r = 0.04, dyn_er = 0),
  delta = 1 / 32,
  start = list(
    formula(y ~ 0.3)
  ),
  state_eq = list(
    formula(y ~ r * y * (1 - (y / k)))
  ),
  meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, 1))),
  dynamic_error = list(
    formula(y ~ dyn_er)
  ),
  boundary = function(x) ifelse(x >= 0, x, 0)
)

sim <- sim_tsm(latent_change)
data <- get_tsm_data(sim, T, T, T)

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

dynr_fit <- tryCatch(dynr.cook(model),
  error = function(cond) {
    dynr_fit <- "Optimizer Error"
    return(dynr_fit)
  }
)

summary(dynr_fit)

test_cv <- c()
for (i in seq_len(nrow(data))) {
  temp_data <- data
  temp_data$y_obs[i] <- NA
  dynr_dat <- dynr.data(
    data.frame(
      id = 1,
      time = temp_data$time,
      y_obs = temp_data$y_obs
    ),
    id = "id", time = "time", observed = "y_obs"
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


  coef(model) <- coef(dynr_fit)
  dynr_fit_temp <- dynr.cook(model, optimization_flag = FALSE)
  test_cv <- c(test_cv, data$y_obs[i] - dynr_fit_temp@eta_smooth_final[1, i])
}

errors <- data$y_obs - dynr_fit@eta_smooth_final[1, ]

A_i_i <- as.vector(dynr_fit@error_cov_smooth_final[1, 1, ]) /
  coef(dynr_fit)["meas_er"]

errors
test_cv

sum(A_i_i - (1 - (errors / test_cv)))

nrow(data) * sum(errors^2) / (nrow(data) - sum((1 - (errors / test_cv))))^2
nrow(data) * sum(errors^2) / (nrow(data) - sum(A_i_i))^2

### Trying to fit the cusp model -----------------------------------------------
# ToDo: Add missing data to improve model fitting?

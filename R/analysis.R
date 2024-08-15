#' ----------------------------------------------------------------------------#
#' Title: Main Analyisis                                                       #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 12-04-2024                                                    #
#' -----                                                                       #
#' Last Modified: 06-08-2024                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#

### Dependencies ---------------------------------------------------------------
library(ggplot2)
library(ggdist)
library(data.table)

# Load functions
invisible(sapply(
  c(paste0("./R/helper/", dir(path = "./R/helper"))),
  source
))

### Load data ------------------------------------------------------------------
load("R/data/simulation_data_17_05_2024_06_59.Rdata")
load("R/data/simulation_results_17_05_2024_06_59.Rdata")
load("R/data/err_test_data_27_05_2024_15_11.Rdata")
load("R/data/err_test_results_27_05_2024_15_11.Rdata")
load("R/data/exp_growth_test_results_28_05_2024_04_13.Rdata")
load("R/data/gam_test_data_06_08_2024_15_56.Rdata")
load("R/data/gam_test_results_06_08_2024_15_56.Rdata")

res <- as.data.table(res)
sim <- as.data.table(sim)

res[, weeks := ifelse(time == 100, 1, 2)]
res[, dyn_var := dyn_er^2]
res[, meas := 1 / (stepsize * (7 / 50))]


res$model[res$model == " exp_growth"] <- "Exponential Growth"
res$model[res$model == "log_growth"] <- "Logistic Growth"
res$model[res$model == "damped_oscillator"] <- "Damped Oscillator"
res$model[res$model == "cusp_catastrophe"] <- "Cusp Catastrophe"

res$method[res$method == "locpol"] <- "Local Polynomial Regression"
res$method[res$method == "gp"] <- "Gaussian Process Regression"
res$method[res$method == "gam"] <- "Generall Additive Modelling"
res$method[res$method == "dynm"] <- "Dynamic Modelling"

res$model <- factor(res$model,
  levels = c(
    "Exponential Growth", "Logistic Growth",
    "Damped Oscillator", "Cusp Catastrophe"
  )
)

res$method <- factor(res$method,
  levels = c(
    "Local Polynomial Regression", "Gaussian Process Regression",
    "Generall Additive Modelling", "Dynamic Modelling"
  )
)

### Analysis -------------------------------------------------------------------
## Descriptives
res_summary <- res[, .(
  mse_mean = mean(mse, na.rm = TRUE),
  mse_se = sd(mse, na.rm = TRUE) / sqrt(.N),
  mse_missing = sum(is.na(mse)),
  gcv_mean = mean(gcv, na.rm = TRUE),
  gcv_se = sd(gcv, na.rm = TRUE) / sqrt(.N),
  gcv_missing = sum(is.na(gcv)),
  ci_coverage_mean = mean(ci_coverage, na.rm = TRUE),
  ci_coverage_se = sd(ci_coverage, na.rm = TRUE) / sqrt(.N),
  ci_coverage_missing = sum(is.na(ci_coverage))
),
by = .(method, model, time, stepsize, dyn_er)
]

View(res_summary)

## Visulization
x11()
plot_results(res = res, "mse", "mean", "weeks", "meas", "dyn_var")
plot_results(res = res, "gcv", "all", "weeks", "meas", "dyn_var")
plot_results(res = res, "ci_coverage", "all", "weeks", "meas", "dyn_var")

## Clean datas
# Remove all oversmoothed GAMS
ind_gam <- sim[, .I[sapply(method, function(x) {
  est <- x[[3]]@estimate
  ind <- seq_along(est)
  fit <- lm(est ~ ind)
  sum(residuals(fit)^2)
}) < 1e-10]] * 4 - 1 # Find all linear gams n = 773

ind_gam <- sim[, .I[sapply(method, function(x) {
  est <- x[[3]]@estimate
  ind <- seq_along(est)
  fit <- lm(est ~ ind)
  summary(fit)$sigma
}) < 0.14]] * 4 - 1 # Find all gams with insignificant smooth terms n = 847 + 3
ind_gam <- c(ind_gam, c(3028, 6772, 6103, 5521) * 4 - 1) # Remove rest manually

# Remove all undersmoothed GAMS
ind_gam <- c(ind_gam, which(res$gcv == 0)) # Find all gams with gcv = 0

# Remove all dynm with NA confidence interval
ind_dynm <- sim[, .I[sapply(method, function(x) {
  any(sapply(x[[4]]@ci, anyNA))
})]] * 4

# Remove all dynm with weird GCVs
ind_dynm <- c(ind_dynm, na.omit(res[
  method == "dynm",
  .I[gcv > 10 | gcv < .5]
] * 4))

res_clean <- res
res_clean[c(ind_gam, ind_dynm), c("mse", "gcv", "ci_coverage")] <- NA

plot_results(res = res_clean, "mse", "mean", "meas")

dat_id <- res_clean[
  method == "dynm",
  order(gcv, decreasing = TRUE)[1:137]
]
inspect(dat_id[1], "dynm", sim)

data <- sim$dat[[dat_id[20, V1]]]
fit <- gam(y_obs ~ s(time, bs = "tp", k = nrow(data) - 1), data = data)
summary(fit)
plot(fit)

x <- sim$method[[ind]][[3]]@estimate
y <- seq_along(x)
test_fit <- lm(x ~ y)

## Inference
predictor <- c("method", "model", "time", "stepsize", "dyn_er")

# Model selection - Univariate
# mse
fit_fun_mse <- function(rhs, data) {
  form <- as.formula(paste0("mse ~ ", rhs))
  fit <- lm(form, data = data)
  return(fit)
}

mse_sel <- model_search(fit_fun_mse, predictor, data = res_clean)
plot(seq_len(nrow(mse_sel)), mse_sel$w_AIC, type = "l")
plot(seq_len(nrow(mse_sel)), mse_sel$w_BIC, type = "l")

# GCV
fit_fun_gcv <- function(rhs, data) {
  form <- as.formula(paste0("gcv ~ ", rhs))
  fit <- lm(form, data = data)
  return(fit)
}

gcv_sel <- model_search(fit_fun_gcv, predictor, data = res_clean)
plot(seq_len(nrow(gcv_sel)), gcv_sel$w_AIC, type = "l")
plot(seq_len(nrow(gcv_sel)), gcv_sel$w_BIC, type = "l")

# Model fitting
# AIC
fit_aic_mse <- lm(
  as.formula(
    paste0("mse ~ ", mse_sel$predictors[which.max(mse_sel$w_AIC)])
  ),
  data = res_clean
)
anova(fit_aic_mse)

emmeans(fit_aic_mse, specs = c("model", "method", "time"))

fit_aic_gcv <- lm(
  as.formula(
    paste0("gcv ~ ", gcv_sel$predictors[which.max(gcv_sel$w_AIC)])
  ),
  data = res
)
anova(fit_aic_gcv)

# BIC
fit_bic_mse <- lm(
  as.formula(
    paste0("mse ~ ", mse_sel$predictors[which.max(mse_sel$w_BIC)])
  ),
  data = res
)
anova(fit_bic_mse)

fit_bic_gcv <- lm(
  as.formula(
    paste0("gcv ~ ", gcv_sel$predictors[which.max(gcv_sel$w_BIC)])
  ),
  data = res
)
anova(fit_bic_gcv)

### Trying to fit the cusp model -----------------------------------------------
# ToDo: Add missing data to improve model fitting?

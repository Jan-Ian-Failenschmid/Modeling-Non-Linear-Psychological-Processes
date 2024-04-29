#' ----------------------------------------------------------------------------#
#' Title: Main Analyisis                                                       #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 12-04-2024                                                    #
#' -----                                                                       #
#' Last Modified: 29-04-2024                                                   #
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
load("R/data/simulation_results_29_04_2024_08_26.Rdata")

### Analysis -------------------------------------------------------------------
## Descriptives
res_summary <- res %>%
  group_by(method, model, time, stepsize, dyn_er) %>%
  summarise(
    mean_mse = mean(mse, na.rm = TRUE),
    missing_mse = sum(is.na(mse)),
    sd_mse = sd(mse, na.rm = TRUE),
    mean_gcv = mean(gcv, na.rm = TRUE),
    sd_gcv = sd(gcv, na.rm = TRUE),
    missing_gcv = sum(is.na(gcv)),
    mean_ci_coverage = mean(ci_coverage, na.rm = TRUE),
    sd_ci_coverage = sd(ci_coverage, na.rm = TRUE),
    missing_ci_coverage = sum(is.na(ci_coverage))
  )

View(res_summary)

## Visulization
x11()
plot_results(res = res, "mse", time, stepsize)
plot_results(res = res, "gcv", time, stepsize)
plot_results(res = res, "ci_coverage", time, stepsize)

## Monte carlo Error
nsim <- seq(30, 1000, 5)
mc_sd_max <- res_summary %>%
  ungroup() %>%
  select(starts_with("sd")) %>%
  sapply(max, na.rm = TRUE)
mcse <- lapply(
  mc_sd_max, function(mce, nsim) mce / sqrt(nsim),
  nsim = nsim
)
plot(0, 0,
  xlim = c(min(nsim), max(nsim)),
  ylim = c(0, max(sapply(mcse, max))),
  main = "Expected Monte Carlo Error by Simulation Replication",
  xlab = "Simulation Replications",
  ylab = "Minimum Monte Carlo"
)
mapply(
  function(mcse, col, nsim) {
    lines(nsim, mcse, col = col)
  },
  mcse = mcse, col = c("red", "blue", "green"),
  MoreArgs = list(nsim = nsim)
)

legend(700, 0.15,
  legend = c("GCV", "CI-Coverage", "MSE"),
  col = c("blue", "green", "red"), lty = c(1, 1, 1)
)

abline(h = 0.03, lty = 2)

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
    paste0("mse ~ ", mse_sel$predictors[which.max(mse_sel$w_AIC)])
  ),
  data = res
)
anova(fit_aic_mse)

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

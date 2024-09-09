#' ----------------------------------------------------------------------------#
#' Title: Main Analyisis                                                       #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 12-04-2024                                                    #
#' -----                                                                       #
#' Last Modified: 09-09-2024                                                   #
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
library(lmtest)
library(car)
library(effectsize)
library(marginaleffects)
library(emmeans)
library(patchwork)

# Load functions
invisible(sapply(
  c(paste0("./R/helper/", dir(path = "./R/helper"))),
  source
))

### Load data ------------------------------------------------------------------
load("R/data/simulation_data_17_05_2024_06_59.Rdata")
load("R/data/simulation_results_17_05_2024_06_59.Rdata")
load("R/data/err_test_data_27_05_2024_15_11.Rdata")
load("R/data/pilot_data_09_09_2024_02_54.Rdata")
load("R/data/pilot_results_09_09_2024_02_54.Rdata")
load("R/data/pilot_data_09_09_2024_01_47.Rdata")
load("R/data/pilot_results_09_09_2024_01_47.Rdata")

res <- as.data.table(res)
sim <- as.data.table(sim)

res[, weeks := ifelse(time == 100, 1, 2)]
res[, dyn_var := dyn_er^2]
res[, meas := 1 / (stepsize * (7 / 50))]


res$model[res$model == "exp_growth"] <- "Exponential Growth"
res$model[res$model == "log_growth"] <- "Logistic Growth"
res$model[res$model == "damped_oscillator"] <- "Damped Oscillator"
res$model[res$model == "cusp_catastrophe"] <- "Cusp Catastrophe"

res$method[res$method == "locpol"] <- "Local Polynomial Regression"
res$method[res$method == "gp"] <- "Gaussian Process Regression"
res$method[res$method == "gam"] <- "General Additive Modelling"
res$method[res$method == "dynm"] <- "Dynamic Modelling"
res$method[res$method == "simple"] <- "Linear Regression"
res$method[res$method == "poly"] <- "Polynomial Regression"

res$model <- factor(res$model,
  levels = c(
    "Exponential Growth", "Logistic Growth",
    "Damped Oscillator", "Cusp Catastrophe"
  )
)
contrasts(res$model) <- contr.sum(levels(res$model))

res$method <- factor(res$method,
  levels = c(
    "Local Polynomial Regression", "Gaussian Process Regression",
    "General Additive Modelling", "Dynamic Modelling",
    "Linear Regression", "Polynomial Regression"
  )
)
contrasts(res$method) <- contr.sum(levels(res$method))

res$weeks <- factor(res$weeks)
contrasts(res$weeks) <- contr.sum(levels(res$weeks))

res$meas <- factor(res$meas)
contrasts(res$meas) <- contr.sum(levels(res$meas))

res$dyn_var <- factor(res$dyn_var,
  levels = c(.5, 1, 2),
  labels = c(".5", "1", "2")
)
contrasts(res$dyn_var) <- contr.sum(levels(res$dyn_var))

res <- res[method != "Linear Regression", ]

### Descriptives ---------------------------------------------------------------
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

### Visulization ---------------------------------------------------------------
dpi <- 300
x11()
for (what in c("all", "missing")) {
  p1 <- plot_results(res = res, "mse", what, "weeks", "meas", "dyn_var")
  p1 <- p1 + theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
    labs(y = "MSE", fill = "Process", color = "Process")

  p2 <- plot_results(res = res, "gcv", what, "weeks", "meas", "dyn_var")
  p2 <- p2 + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
    labs(y = "GCV", fill = "Process", color = "Process")

  p3 <- plot_results(res = res, "ci_coverage", what, "weeks", "meas", "dyn_var")
  p3 <- p3 + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 5)
  ) +
    labs(
      x = "Simulation Conditions", y = "CI Coverage", ,
      fill = "Process", color = "Process"
    ) +
    ylim(0, 1) +
    annotate("text", x = -1, y = -3, label = "Test") +
    coord_cartesian(clip = "off")

  p_comb <- p1 / p2 / p3 +
    plot_layout(guides = "collect") & theme(legend.position = "right")
  p_comb

  ggsave(paste0("figures/complete_results_", what, ".png"), p_comb,
    width = 1920 * (dpi / 72), height = 1080 * (dpi / 72), units = "px", dpi = dpi, limitsize = FALSE
  )
}

### Clean data -----------------------------------------------------------------
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

### Plot mean results ----------------------------------------------------------
for (var in c(1, 2, 3)) {
  outcome <- switch(var,
    a = "mse",
    b = "gcv",
    c = "ci_coverage"
  )
  var_name <- switch(var,
    a = "MSE",
    b = "GCV",
    c = "CI coverage"
  )
  y_lab_str <- paste0("Mean ", var_name)
  p1 <- plot_results(res = res, outcome, "mean")
  p1 <- p1 + theme(
    axis.text.x = element_blank()
  ) + labs(
    title = "a) Effect of Analysis Method for each Process",
    y = y_lab_str, fill = "Process", color = "Process",
    x = "Process"
  )

  p2 <- plot_results(res = res, outcome, "mean", "weeks")
  p2 <- p2 + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + labs(
    title = "b) Effect of Measurement Period for each Method and Process",
    y = y_lab_str, fill = "Process", color = "Process",
    x = "Measurement Period (in weeks)"
  )

  p3 <- plot_results(res = res, outcome, "mean", "meas")
  p3 <- p3 + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
  ) + labs(
    title = "c) Effect of Measurement Frequency for each Method and Process",
    y = y_lab_str, fill = "Process", color = "Process",
    x = "Measurement Frequency (per day)"
  )

  p4 <- plot_results(res = res, outcome, "mean", "dyn_var")
  p4 <- p4 + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + labs(
    title =
      "d) Effect of Dynamic Error Variance for each Method and Process",
    y = y_lab_str, fill = "Process", color = "Process",
    x = "Dynamic Error Variance"
  )

  p_comb <- p1 / p2 / p3 / p4 +
    plot_layout(guides = "collect") & theme(legend.position = "right")

  p_ranges_y <- c(
    ggplot_build(p_comb[[1]])$layout$panel_scales_y[[1]]$range$range,
    ggplot_build(p_comb[[2]])$layout$panel_scales_y[[1]]$range$range
  )

  if (var == 3) p_ranges_y <- c(0, 1)

  p_comb <- p_comb & ylim(min(p_ranges_y), max(p_ranges_y))

  if (var == 3) {
    p_comb <- p_comb & annotate("rect",
      xmin = -Inf, xmax = Inf, ymin = 0.93, ymax = 0.97,
      alpha = 0.5
    ) & geom_hline(yintercept = 0.95, color = "black")
  }

  ggsave(paste0("figures/mean_results_", outcome, ".png"), p_comb,
    width = 1920 * (dpi / 72), height = 1080 * (dpi / 72), units = "px",
    dpi = dpi, limitsize = FALSE
  )
}

# Mischellaneous
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
predictor <- c("method", "model", "weeks", "meas", "dyn_var")

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
qqnorm(resid(fit_aic_mse))
qqline(resid(fit_aic_mse))
hist(resid(fit_aic_mse))
bptest(fit_aic_mse)

aov_mse <- Anova(fit_aic_mse, singular.ok = TRUE)
effectsize(aov_mse, partial = TRUE)

emmeans(fit_aic_mse, specs = "model")
emmeans(fit_aic_mse, specs = "method")

fit_aic_gcv <- lm(
  as.formula(
    paste0("gcv ~ ", gcv_sel$predictors[which.max(gcv_sel$w_AIC)])
  ),
  data = res
)

qqnorm(resid(fit_aic_gcv))
qqline(resid(fit_aic_gcv))
hist(resid(fit_aic_gcv))
bptest(fit_aic_gcv)

aov_gcv <- Anova(fit_aic_gcv, singular.ok = TRUE)
effectsize(aov_gcv, partial = TRUE)
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

ilustr <- sim[, model_name := sapply(gen_model, function(x) {
  x@model_name
})][,
  .I[time == 200 & stepsize == unique(stepsize)[1] & dyn_er == sqrt(1)],
  by = model_name
][, .SD[2], by = model_name]

method <-
  c(
    "Local Polynomial Regression",
    "Gaussian Process Regression",
    "General Additive Model",
    "Simple Linear Regression",
    "Polynomial Regression"
  )

process <- c(
  "Exponential Growth", "Logistic Growth",
  "Damped Oscillator", "Cusp Catastrophe"
)
png(
  file = "./figures/smooth.png",
  width = 1960, height = 1080
)
par(mfrow = c(4:5))
for (j in 1:4) {
  for (i in 1:5) {
    plot(sim$method[[ilustr$V1[j]]][[i]],
      sim = sim,
      axes = FALSE, xlab = "", ylab = ""
    )
    axis(1,
      at = c(0, 50, 100, 150, 200),
      labels = c("", "Week 1", "", "Week 2", ""), line = c(0, 200)
    )
    if (j == 1) {
      title(method[i])
    } else if (j == 4) {
      title(xlab = "Time", main = NULL)
    }
    if (i == 1) {
      title(ylab = process[i], main = NULL)
    } else {
      title(main = NULL)
    }
  }
}
dev.off()

#' ----------------------------------------------------------------------------#
#' Title: Main Analyisis                                                       #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 12-04-2024                                                    #
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

### Dependencies ---------------------------------------------------------------
library(ggplot2)
library(ggdist)
library(data.table)
library(papaja)
library(lmtest)
library(car)
library(effectsize)
# library(marginaleffects)
# library(emmeans)
library(patchwork)

# Load functions
invisible(sapply(
  c(paste0("./R/helper/", dir(path = "./R/helper"))),
  source
))

### Pilot data -----------------------------------------------------------------
load("R/data/old_data/pilot_results_09_09_2024_01_47.Rdata")
res1 <- as.data.table(res) # LPR, GP, GAM
res1 <- res1[!method %in% c("dynm", "simple", "poly"), ]

load("R/data/old_data/pilot_results_10_09_2024_14_31.Rdata")
res2 <- as.data.table(res) # Parametric models
res2 <- res2[!method %in% c("gam"), ]

load("R/data/old_data/pilot_results_11_09_2024_18_39.Rdata")
res3 <- as.data.table(res) # Parametric models
res3 <- res3[!method %in% c("gam"), ]

res_pilot <- rbind(res1, res2, res3)

res_summary <- res_pilot[, .(
  mse_mean = mean(mse, na.rm = TRUE),
  mse_sd = sd(mse, na.rm = TRUE),
  mse_missing = sum(is.na(mse)),
  gcv_mean = mean(gcv, na.rm = TRUE),
  gcv_sd = sd(gcv, na.rm = TRUE),
  gcv_missing = sum(is.na(gcv)),
  ci_coverage_mean = mean(ci_coverage, na.rm = TRUE),
  ci_coverage_sd = sd(ci_coverage, na.rm = TRUE),
  ci_coverage_missing = sum(is.na(ci_coverage))
),
by = .(method, model, time, stepsize, dyn_er)
]

nsim <- seq(30, 1000, 5)

# Extract the largest standard deviation by metric across conditions
mc_sd_max <- sapply(res_summary[, .(mse_sd, gcv_sd, ci_coverage_sd)],
  max,
  na.rm = TRUE
)

# Devide the largest stardard deviations by the sample sizes to find the
# MC standard errors
mcse <- lapply(
  mc_sd_max, function(mce, nsim) mce / sqrt(nsim),
  nsim = nsim
)

### Load and prepare data ------------------------------------------------------
load("R/data/combined_results_06_02_2025_21_17.Rdata")
res <- as.data.table(res)
res <- res[!method %in% c("simple", "poly_orth"), ]


# load("R/data/old_data/pilot_data_09_09_2024_01_47.Rdata")
# sim1 <- as.data.table(sim) # LPR, GP, GAM

# load("R/data/old_data/pilot_data_10_09_2024_14_31.Rdata")
# sim2 <- as.data.table(sim) # Parametric models

# load("R/data/old_data/pilot_data_11_09_2024_18_39.Rdata")
# sim3 <- sim # Correlated Polynomial Regression

load("R/data/combined_data_06_02_2025_21_17.Rdata")

res[, SP := ifelse(time == 100, 0.5, 1)]
res[, DEV := dyn_er^2]
res[, SF := (1 / (stepsize * (7 / 50))) / 3]


res$model[res$model == "exp_growth"] <- "Exponential Growth"
res$model[res$model == "log_growth"] <- "Logistic Growth"
res$model[res$model == "damped_oscillator"] <- "Damped Oscillator"
res$model[res$model == "cusp_catastrophe"] <- "Cusp Catastrophe"

res$method[res$method == "locpol"] <- "Local Polynomial Regression"
res$method[res$method == "gp"] <- "Gaussian Process Regression"
res$method[res$method == "gam"] <- "Generalized Additive Modelling"
res$method[res$method == "dynm"] <- "Parametric Modelling"
res$method[res$method == "simple"] <- "Linear Regression"
res$method[res$method == "poly"] <- "Global Polynomial Regression"

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
    "Generalized Additive Modelling",
    "Linear Regression", "Global Polynomial Regression", "Parametric Modelling"
  )
)
contrasts(res$method) <- contr.sum(levels(res$method))

res$SP <- factor(res$SP)
contrasts(res$SP) <- contr.sum(levels(res$SP))

res$SF <- factor(res$SF)
contrasts(res$SF) <- contr.sum(levels(res$SF))

res$DEV <- factor(res$DEV,
  levels = c(.5, 1, 2),
  labels = c(".5", "1", "2")
)
contrasts(res$DEV) <- contr.sum(levels(res$DEV))

### Descriptives ---------------------------------------------------------------
res_summary <- res[, .(
  mse_mean = mean(mse, na.rm = TRUE),
  mse_se = sd(mse, na.rm = TRUE) / sqrt(.N),
  mse_missing = sum(is.na(mse)) / .N,
  gcv_mean = mean(gcv, na.rm = TRUE),
  gcv_se = sd(gcv, na.rm = TRUE) / sqrt(.N),
  gcv_missing = sum(is.na(gcv)),
  ci_coverage_mean = mean(ci_coverage, na.rm = TRUE),
  ci_coverage_se = sd(ci_coverage, na.rm = TRUE) / sqrt(.N),
  ci_coverage_missing = sum(is.na(ci_coverage))
),
by = .(method)
]

# View(res_summary)

### Visulization ---------------------------------------------------------------
dpi <- 300

# Complete results plot
p1 <- plot_results(res = res, "mse", "all", "SP", "SF", "DEV")
p1 <- p1 +
  theme_apa() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(y = "MSE", fill = "Process", color = "Process")

p2 <- plot_results(res = res, "gcv", "all", "SP", "SF", "DEV")
p2 <- p2 +
  theme_apa() +
  theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(y = "GCV", fill = "Process", color = "Process")

p3 <- plot_results(res = res, "ci_coverage", "all", "SP", "SF", "DEV")
p3 <- p3 +
  theme_apa() +
  theme(
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 5)
  ) +
  labs(
    x = "Simulation Conditions", y = "CI Coverage", ,
    fill = "Process", color = "Process"
  ) +
  ylim(0, 1)

p_comb <- p1 / p2 / p3 +
  plot_layout(guides = "collect") & theme(legend.position = "right")
p_comb

ggsave("figures/complete_results_all.png", p_comb,
  width = 1920 * (dpi / 72), height = 1080 * (dpi / 72), units = "px",
  dpi = dpi, limitsize = FALSE
)

# Missing data plot
pmiss <- plot_results(res = res, "mse", "missing", "SP", "SF", "DEV")
pmiss <- pmiss +
  scale_y_continuous(expand = expansion(mult = 0, add = 0)) +
  theme_apa() +
  theme(
    axis.text.x = element_text(size = 5)
  ) +
  labs(y = "Proportion Not-converged", fill = "Process", color = "Process")

ggsave("figures/complete_results_missing.png", pmiss,
  width = 1920 * (dpi / 72), height = 1080 * (dpi / 72), units = "px", dpi = dpi, limitsize = FALSE
)

### Plot mean results ----------------------------------------------------------
# Simulation conditions plots
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
  p1 <- p1 +
    theme_apa() +
    theme(
      axis.text.x = element_blank(),
      text = element_text(size = 26)
    ) + labs(
      title = "a) Effect of Analysis Method for each Process",
      y = y_lab_str, fill = "Process", color = "Process",
      x = "Process"
    )

  p2 <- plot_results(res = res, outcome, "mean", "SP")
  p2 <- p2 +
    theme_apa() +
    theme(
      strip.text.x = element_blank(),
      text = element_text(size = 26)
    ) + labs(
      title = "c) Effect of Measurement Period for each Method and Process",
      y = y_lab_str, fill = "Process", color = "Process",
      x = "Measurement Period"
    )

  p3 <- plot_results(res = res, outcome, "mean", "SF")
  p3 <- p3 +
    theme_apa() +
    theme(
      strip.text.x = element_blank(),
      text = element_text(size = 26)
    ) + labs(
      title = "d) Effect of Measurement Frequency for each Method and Process",
      y = y_lab_str, fill = "Process", color = "Process",
      x = "Measurement Frequency"
    )

  p4 <- plot_results(res = res, outcome, "mean", "DEV")
  p4 <- p4 +
    theme_apa() +
    theme(
      strip.text.x = element_blank(),
      text = element_text(size = 26)
    ) + labs(
      title =
        "b) Effect of Dynamic Error Variance for each Method and Process",
      y = y_lab_str, fill = "Process", color = "Process",
      x = "Dynamic Error Variance"
    )

  p_comb <- p1 / p4 / p2 / p3 +
    plot_layout(guides = "collect") & theme(legend.position = "right")

  p_ranges_y <- c(
    ggplot_build(p_comb[[1]])$layout$panel_scales_y[[1]]$range$range,
    ggplot_build(p_comb[[2]])$layout$panel_scales_y[[1]]$range$range,
    ggplot_build(p_comb[[3]])$layout$panel_scales_y[[1]]$range$range,
    ggplot_build(p_comb[[4]])$layout$panel_scales_y[[1]]$range$range
  )

  if (var == 3) p_ranges_y <- c(0, 1)

  p_comb <- p_comb & ylim(min(p_ranges_y), max(p_ranges_y))

  if (var == 3) {
    p_comb <- p_comb & annotate("rect",
      xmin = -Inf, xmax = Inf, ymin = 0.89, ymax = 1,
      alpha = 0.3
    )
  }

  ggsave(paste0("figures/mean_results_", outcome, ".png"), p_comb,
    width = 1920 * (dpi / 72), height = 1080 * (dpi / 72), units = "px",
    dpi = dpi, limitsize = FALSE
  )
}


# Smoothing plot
ilustr <- sim[, model_name := sapply(gen_model, function(x) {
  x@model_name
})][,
  .I[time == 200 & stepsize == unique(stepsize)[1] & dyn_er == sqrt(1)],
  by = model_name
][, .SD[3], by = model_name]

method <-
  c(
    "Local Polynomial Regression",
    "Gaussian Process Regression",
    "Generalized Additive Model",
    "Parametric Modelling",
    "",
    "Global Polynomial Regression"
  )

process <- c(
  "Exponential Growth", "Logistic Growth",
  "Damped Oscillator", "Cusp Catastrophe"
)
png(
  file = "./figures/smooth.png",
  width = 1960, height = 1080
)

par(
  mfrow = c(4, 5), cex.lab = 2.5, cex.main = 2.5, cex.axis = 2.5,
  mar = c(6, 6, 4, 4) + .1, lwd = 2.5
)

for (j in 1:4) {
  for (i in c(1:3, 6, 4)) {
    plot(sim$method[[ilustr$V1[j]]][[i]],
      sim = sim,
      axes = FALSE, xlab = "", ylab = ""
    )

    axis(2, at = seq(-10, 10, 2), cex.axis = 2.5)

    axis(1,
      at = c(0, 50, 100, 150, 200),
      labels = c("", "First Half", "", "Sec. Half", ""), cex.axis = 2.5,
      padj = 1
    )

    if (j == 1) {
      title(method[i])
    } else if (j == 4) {
      title(xlab = "Time", main = NULL, line = 5)
    }
    if (i == 1) {
      title(ylab = process[j], main = NULL)
    } else {
      title(main = NULL)
    }
  }
}
dev.off()

### ANOVA ----------------------------------------------------------------------
# MSE
res_aov <- res[method != "Parametric Modelling", c(1:4, 9:11)]
res_aov$method <- factor(res_aov$method,
  levels = c(
    "Local Polynomial Regression", "Gaussian Process Regression",
    "Generalized Additive Modelling", "Global Polynomial Regression"
  )
)
contrasts(res_aov$method) <- contr.sum(levels(res_aov$method))

names(res_aov) <- c("Method", "Process", "mse", "gcv", "SP", "DEV", "SF")

# MSE
mse_lm <- lm(mse ~ Method * Process * SP * SF * DEV, data = res_aov)

qqnorm(resid(mse_lm))
qqline(resid(mse_lm))
hist(resid(mse_lm), probability = TRUE)
curve(dnorm(x, sd = sd(resid(mse_lm))), add = TRUE, yaxt = "n")
psych::describe(resid(mse_lm))

# shapiro.test(df$Gewicht)
bptest(mse_lm)
leveneTest(mse_lm)

mse_aov <- Anova(mse_lm, type = "III")
mse_aov_test <- Anova(mse_lm, type = "III", white.adjust = "hc3")


apa_table(mse_aov_test)

mse_etas <- eta_squared(mse_aov)

# GCV
gcv_lm <- lm(gcv ~ Method * Process * SP * SF * DEV, data = res_aov)

qqnorm(resid(gcv_lm))
qqline(resid(gcv_lm))
hist(resid(gcv_lm), probability = TRUE)
curve(dnorm(x, sd = sd(resid(gcv_lm))), add = TRUE, yaxt = "n")
psych::describe(resid(gcv_lm))

# shapiro.test(df$Gewicht)
bptest(gcv_lm)
leveneTest(gcv_lm)

gcv_aov <- Anova(gcv_lm, type = "III")
gcv_aov_test <- Anova(gcv_lm, type = "III", white.adjust = "hc3")

gcv_etas <- eta_squared(gcv_aov)

# Table outputs
etsq_tabel <- as.data.table(cbind(mse_etas[, 1:2], gcv_etas[, 2]))
names(etsq_tabel) <- c("Effect", "MSE Eta2_partial", "GCV Eta2_partial")
ind <- as.vector(etsq_tabel[, 2] >= 0.01 | etsq_tabel[, 3] >= 0.01)

writeLines(
  apa_table(etsq_tabel[ind, ]),
  "./tables/eta_squares.txt"
)

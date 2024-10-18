#' ----------------------------------------------------------------------------#
#' Title: Real Data Demonstration                                              #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 23-05-2024                                                    #
#' -----                                                                       #
#' Last Modified: 18-10-2024                                                   #
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
# install.packages("ggdist") # Plotting
# install.packages("data.table") # Data table for storing the simulation grid
# install.packages("ggplot2") # Plotting
# install.packages("patchwork")
# install.packages("papaja") # Export package list

library(mgcv) # GAM's
library(ggdist) # Plotting
library(data.table) # Data table for storing the simulation grid
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

# Time since fist assessment
df_raw[, time0 := difftime(time, min(time), units = "hours"), by = UUID]

df_raw$UUID <- as.factor(df_raw$UUID)
levels(df_raw$UUID) <- 1:118
df <- df_raw[UUID != unique(UUID)[110], ] # Part 110 excluded due to no var

# Group data by UUID
grouped_df <- df[, .(data = list(.SD)),
  by = UUID, .SDcols = c("DEP_ES", "time", "time0")
]

# Descriptives
nrow(grouped_df) # Number of participants
summary(df[, .(age = unique(AGE_BL)), by = UUID][, age])
summary(df[, .(gender = unique(as.factor(GENDER_BL))), by = UUID][, gender])


### Fit GAMs -------------------------------------------------------------------
# Fit gams
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

# Extract predictions
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

# Extract edf
grouped_df[, gam_edf := lapply(gam_fit, function(fit) {
  summary(fit)$edf
})]

# Extract residuals
grouped_df[, gam_resid := mapply(
  function(data, gam_pred) {
    na.omit(data$DEP_ES) - gam_pred$est
  },
  data = data, gam_pred = gam_pred, SIMPLIFY = FALSE
)]

# Calculate mse
grouped_df[, gam_mse := lapply(
  gam_resid,
  function(resid) mean(resid^2, na.rm = TRUE)
)]

# Calculate gcv
grouped_df[, gam_gcv := mapply(
  function(fit, resid) {
    n <- length(resid)
    n * sum(resid^2) / (n - sum(influence(fit)))^2
  },
  fit = gam_fit, resid = gam_resid
)]

# Calculate simga
grouped_df[, gam_sigma := lapply(
  gam_fit,
  function(fit) {
    sqrt(sum(fit$residuals^2) / fit$df.residual)
  }
)]

# Make high resolution predictions for plotting
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
# Gam with the lowest edf
np <- 10
sum(grouped_df$gam_edf < 1.001)
gam_low <- plot_ml_pred(
  id = grouped_df$UUID[order(unlist(grouped_df$gam_edf))[1:np]],
  pred = grouped_df$gam_pred_plot[order(unlist(grouped_df$gam_edf))[1:np]],
  dat = grouped_df$data[order(unlist(grouped_df$gam_edf))[1:np]],
  var = "DEP_ES", time = "time0"
)
gam_low <- gam_low +
  ggtitle("a) Ten Participants with the lowest EDF") +
  xlab("") +
  ylab("Momentary depression")

# Gam with the highest edf
gam_high <- plot_ml_pred(
  id = grouped_df$UUID[order(unlist(grouped_df$gam_edf),
    decreasing = TRUE
  )[1:np]],
  pred = grouped_df$gam_pred_plot[order(unlist(grouped_df$gam_edf),
    decreasing = TRUE
  )[1:np]],
  dat = grouped_df$data[order(unlist(grouped_df$gam_edf),
    decreasing = TRUE
  )[1:np]],
  var = "DEP_ES", time = "time0"
)
gam_high <- gam_high +
  ggtitle("b) Ten Participants with the highest EDF") +
  xlab("") +
  ylab("")

# unique(ggplot_build(gam_high)$data[[1]]$fill)

# Oscillating GAM 1
gam_osc1 <- plot_ml_pred(
  id = grouped_df$UUID[12],
  pred = grouped_df$gam_pred_plot[12],
  dat = grouped_df$data[12],
  var = "DEP_ES", time = "time0"
)
gam_osc1 <- gam_osc1 +
  ggtitle("c) Smooth for Participant 61") +
  xlab("Time since first assessment") +
  ylab("Momentary depression") +
  scale_fill_manual(values = "#00BFC4") +
  scale_color_manual(values = "#00BFC4")

# Oscillating GAM 2
gam_osc2 <- plot_ml_pred(
  id = grouped_df$UUID[8],
  pred = grouped_df$gam_pred_plot[8],
  dat = grouped_df$data[8],
  var = "DEP_ES", time = "time0"
)
gam_osc2 <- gam_osc2 +
  ggtitle("d) Smooth for Participant 83") +
  xlab("Time since first assessment") +
  ylab("") +
  scale_fill_manual(values = "#F8766D") +
  scale_color_manual(values = "#F8766D")

complete_plot <- (gam_low | gam_high) /
  (gam_osc1 | gam_osc2)

# Save plot
dpi <- 300
ggsave("figures/demonstration_smooths.png", complete_plot,
  width = 1920 * (dpi / 72), height = 1080 * (dpi / 72), units = "px",
  dpi = dpi, limitsize = FALSE
)

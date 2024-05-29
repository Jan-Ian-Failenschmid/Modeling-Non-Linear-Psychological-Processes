#' ----------------------------------------------------------------------------#
#' Title:                                                                      #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 28-05-2024                                                    #
#' -----                                                                       #
#' Last Modified: 28-05-2024                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#
set.seed(12345)
cex <- 1.5
# Exemplar plots -----

png(
  file = paste0(fig_path, "/damp_osc.png"),
  width = 1960, height = 1080
)

make_exemplar_plot(damp_osc,
  main = NULL,
  xlab = "Time", ylab = "", cex = cex, xaxt = "n",
)

dev.off()

# Non-linearity demonstration ------
plot_pred <- function(data, pred, ind, var, ...) {
  time <- data[[ind]][, time]
  eval_time <- seq(min(time), max(time), length.out = eval_points)
  obs <- unlist((data[[ind]][, var, with = FALSE]))
  plot(x = time, y = obs, ...)
  lines(x = eval_time, y = pred[[ind]]$est, col = "#5D4641")
  lines(x = eval_time, y = pred[[ind]]$ub, lty = 2, col = "#5D4641")
  lines(x = eval_time, y = pred[[ind]]$lb, lty = 2, col = "#5D4641")
}

eval_points <- 200

grouped_df[, gam_fit := lapply(data, function(data) {
  gam(DEP_ES ~ s(as.numeric(time),
    bs = "tp", k = 34
  ), data = data)
})]

grouped_df[, gam_pred := mapply(function(fit, data) {
  eval_time <- data.frame(
    time = seq(min(data[, time]), max(data[, time]),
      length.out = eval_points
    )
  )
  inference <- predict(fit, newdata = eval_time, se.fit = TRUE)
  list(
    est = as.vector(inference$fit),
    ub = as.vector(inference$fit + (qnorm(0.975) * inference$se.fit)),
    lb = as.vector(inference$fit - (qnorm(0.975) * inference$se.fit))
  )
}, fit = gam_fit, data = data, SIMPLIFY = FALSE)]

png(
  file = paste0(fig_path, "/demonstration.png"),
  width = 1960, height = 1080
)

par(mfrow = c(2, 2))
cex <- 1.4
for (i in c(1, 2, 3, 4)) {
  plot_pred(grouped_df[, data], grouped_df[, gam_pred], i, "DEP_ES",
    xlab = "Time", ylab = "Depression", main = paste("Participant", i),
    cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex
  )
}

dev.off()

# Method illustration --------------------
# Locpols
plot_local_polynomial <- function(
    x, y, x_eval, h, p = 1,
    K = dnorm, threshold, alpha = 1) {
  # Create X matrix of the local polinomial smoother (x_i - x)^p
  X <- as.matrix(sapply(0:p, function(p, x, x_eval) (x - x_eval)^p,
    x = x, x_eval = x_eval
  ))

  # Calculate weights for the distances between each x points and x x.eval
  W <- diag(1, length(x)) * (K((x_eval - x) / h))

  # Solve local polynomial equation to get to influence vector
  beta <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y

  # Make predictions for the local polynomial
  pred <- X %*% beta

  lines(x[which(diag(W) > threshold)], pred[which(diag(W) > threshold)],
    col = rgb(
      red = 93, green = 70, blue = 65, alpha = alpha,
      maxColorValue = 255
    )
  )

  points(x[x_eval], pred[x_eval],
    pch = 19,
    col = rgb(
      red = 93, green = 70, blue = 65, alpha = alpha,
      maxColorValue = 255
    ),
    cex = 2
  )

  return(NULL)
}

png(
  file = paste0(fig_path, "/locpol_demonstration7.png"),
  width = 1960, height = 1080
)

plot(dat$time, dat$y_obs,
  main = NULL,
  xlab = "Time", ylab = "Y", xaxt = "n",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex
) # Data points

lines(dat$time, dat$y, lty = 2, lwd = 1) # State function

lines(dat$time, lp_fit$Estimate[, 5], lwd = 1, col = "#5D4641") # Estimated line

points(
  x = c(49, 51, 53, 55), y = lp_fit$Estimate[c(50, 52, 54, 56), 5],
  col = "#5D4641", cex = 2, pch = 19
)

# plot_local_polynomial(
#     x = dat$time, y = dat$y_obs, x_eval = 56,
#     h = lp_fit$Estimate[1, 2], p = lp_fit$opt$p, threshold = 0.02, alpha = 255
# )

dev.off()

# Gaussian Processes
gp_post_pred <- function(
    fit, time, obs, state, f_name, draws = 200,
    alpha = 1, lwd = 1, rect = FALSE,
    quant = c(.60, .80, .95), ...) {
  #' Function for the posterior predictive plotting of a gaussian process

  f_draws <- as.data.frame(fit$draws(f_name, format = "draws_df"))
  f_draws <- f_draws[, seq(1, ncol(f_draws) - 3)]

  if (!rect) {
    plot(x = time, y = obs, ...)

    lapply(sample(seq(1, nrow(f_draws)), draws),
      function(x, f_draws, time, lwd, alpha) {
        lines(
          x = time, y = f_draws[x, seq_len(length(time))],
          col = rgb(
            red = 93, green = 70, blue = 65, alpha = alpha,
            maxColorValue = 255
          ),
          lwd = lwd
        )
      },
      f_draws = f_draws, time = time, alpha = alpha, lwd = lwd
    )

    lines(x = time, y = state)
  } else if (rect) {
    plot_poly <- function(x, quants, time) {
      quant_names <- c((1 - x) / 2, 1 - ((1 - x) / 2))
      quants <- quants[, names(quants) %in% quant_names]
      polygon(
        x = c(time, rev(time)),
        y = c(quants[time, 1], quants[rev(time), 2]),
        col = rgb(red = 1, green = x / 2, blue = 0)
      )
    }

    quants <- lapply(f_draws, function(x, quant) {
      quant <- c((1 - quant) / 2, 1 - ((1 - quant) / 2))
      quant_val <- quantile(x, probs = quant)
      names(quant_val) <- quant
      return(quant_val)
    }, quant = quant)

    quants <- as.data.frame(do.call(rbind, quants))

    plot(x = time, y = obs, type = "n")

    for (x in quant[order(quant, decreasing = TRUE)]) {
      plot_poly(x, quants, time)
    }

    points(x = time, y = obs)

    lines(x = time, y = state)
  }
}

png(
  file = paste0(fig_path, "/gaussian_processes3.png"),
  width = 1960, height = 1080
)

set.seed(42)
gp_post_pred(gp_fit,
  f_name = "f_post_predict", time = dat$time,
  obs = dat$y_obs, state = rep(-20, 201), alpha = 255 / 5, draws = 200,
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex,
  main = NULL, xlab = "Time", ylab = "Y", xaxt = "n"
)
lines(dat$time, dat$y, lty = 2, lwd = 1) # State function

dev.off()

png(
  file = paste0(fig_path, "/gam2.png"),
  width = 1960, height = 1080
)

plot(dat$time, dat$y_obs,
  main = NULL,
  xlab = "Time", ylab = "Y", xaxt = "n",
  cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex
) # Data points

lines(dat$time, dat$y, lty = 2, lwd = 1) # State function

for (i in seq_len(ncol(predmat))) {
  lines(
    x = dat$time, y = predmat[, i],
    col = rgb(
      red = 93, green = 70, blue = 65, alpha = 255,
      maxColorValue = 255
    )
  )
}

# lines(dat$time, fitted(gam_fit), lwd = 1, col = "#5D4641") # Predictions

dev.off()

# Simulation Results
ggsave("ci_results.png",
  device = "png", path = fig_path,
  height = 1080, width = 1960, units = "px", dpi = 100
)

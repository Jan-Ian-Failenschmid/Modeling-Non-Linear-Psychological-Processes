#' ----------------------------------------------------------------------------#
#' Title: Simulation Main File                                                 #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 10-04-2024                                                    #
#' -----                                                                       #
#' Last Modified: 30-08-2024                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#

#' Main simulation file. This file loads in all dependencies, defines the
#' simulation parameters and runs the simulation. The results are automatically
#' extracted and are saved alongside the generated data.
#'
#' This file is ideally executed in a shell by running:
#' Rscript R/simulation.R

### Dependencies ---------------------------------------------------------------
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

# Create data directory
out_dir <- "./R/data"
if (!file.exists(out_dir)) dir.create(out_dir)

# Document Session Info
writeLines(capture.output(sessionInfo()), "./R/sessionInfo.txt")

### Simulation parameters ------------------------------------------------------
## Generative Models
exp_growth <- new("gen_model",
  model_name = "exp_growth",
  model_type = "DE",
  # Set dyn_er to 0 and overwrite during the simulation
  pars = list(yr = 0.02, ya = 2, dyn_er = 0),
  delta = (50 / 7) / 234, # Integer divides 3, 6, and 9 obs. per day
  stepsize = (50 / 7) / 3, # Gets overwritten anyway
  start = list(
    formula(y ~ -2)
  ),
  state_eq = list(
    formula(y ~ yr * ya - yr * y)
  ),
  dynamic_error = list(
    formula(y ~ dyn_er)
  ),
  meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, 1)))
)

log_growth <- new("gen_model",
  model_name = "log_growth",
  model_type = "DE",
  # Set dyn_er to 0 and overwrite during the simulation
  pars = list(k = 4.3, r = 0.04, dyn_er = 0),
  delta = (50 / 7) / 234, # Integer divides 3, 6, and 9 obs. per day
  stepsize = (50 / 7) / 3, # Gets overwritten anyway
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

cusp_catastrophe <- new("gen_model",
  model_name = "cusp_catastrophe",
  model_type = "DE",
  # Set dyn_er to 0 and overwrite during the simulation
  delta = (50 / 7) / 234, # Integer divides 3, 6, and 9 obs. per day
  stepsize = (50 / 7) / 3, # Gets overwritten anyway
  pars = list(a = -5, dyn_er = 0, omega = -(2 * pi / 50)^2),
  start = list(
    formula(y ~ 1.9),
    formula(v ~ 0),
    formula(b ~ -10)
  ),
  state_eq = list(
    formula(y ~ -(4 * y^3 + 2 * a * y + b)),
    formula(b ~ v),
    formula(v ~ omega * b)
  ),
  dynamic_error = list(
    formula(y ~ dyn_er),
    formula(b ~ 0),
    formula(v ~ 0)
  ),
  meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, 1)))
)

damped_oscillator <- new("gen_model",
  model_name = "damped_oscillator",
  model_type = "DE",
  # Set dyn_er to 0 and overwrite during the simulation
  delta = (50 / 7) / 234, # Integer divides 3, 6, and 9 obs. per day
  stepsize = (50 / 7) / 3, # Gets overwritten anyway
  pars = list(k = 0.01, c = 0.1, dyn_er = 0),
  start = list(
    formula(y ~ 2),
    formula(v ~ 0)
  ),
  state_eq = list(
    formula(y ~ v),
    formula(v ~ -2 * k * v - c^2 * y)
  ),
  dynamic_error = list(
    formula(y ~ dyn_er),
    formula(v ~ 0)
  ),
  meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, 1)))
)

## Anayisis Methods
gam <- new("method_gam",
  method_name = "gam"
)

locpol <- new("method_locpol",
  method_name = "locpol"
)

gp <- new("method_gp",
  method_name = "gp"
)

dynm <- new("method_dynm",
  method_name = "dynm"
)

### Run simulation -------------------------------------------------------------
repetitions <- 30 # Number of repetitions in the pilot sample
mc_error_target <- 0.05 # Desired monte carlo error
for (run in c("pilot", "simulation")) {
  # Set seed
  if (run == "pilot") {
    set.seed(12345)
  } else if (run == "simulation") {
    set.seed(54321)
  }

  # Run simulation
  cat("Running", run, "---------------------------------------------------\n\n")
  system.time({
    sim <- simulate(
      gen_model_list = list(
        exp_growth, log_growth, damped_oscillator, cusp_catastrophe
      ),
      method_list = list(locpol, gp, gam, dynm),
      conditions = list(
        time = c(100, 200), # 2 & 4 weeks rescaled to 1 week = 50 units
        # 3, 6, 9, measurements per day
        stepsize = c((50 / 7) / 3, (50 / 7) / 6, (50 / 7) / 9),
        # dynamic error variances of 12.5%, 25%, and 50% of the process range
        dyn_er = sqrt(c(.5, 1, 2))
      ),
      repetitions = repetitions,
      out_dir = out_dir
    )
  })

  # Check object size of sim grid
  cat(utils::object.size(sim)[1] / (1e6), "mb")

  # Save simulation
  sim_time <- format(Sys.time(), "%d_%m_%Y_%H_%M")
  save(sim, file = paste0(out_dir, "/", run, "_data_", sim_time, ".Rdata"))

  # Extract and save results
  res <- extract_results(sim)
  save(res, file = paste0(out_dir, "/", run, "_results_", sim_time, ".Rdata"))

  if (run == "pilot") {
    # rm(sim) # Remove sim from working environment
    res <- as.data.table(res)
    # Calculate summary statistics for each condition
    res_summary <- res[, .(
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

    # Create sequence of sample sizes
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
    # ToDo: Does not currently work
    # Find the largest sample size such that the MC error is smaller than
    # the desired criterion for all three metrics
    n_ind <- max(sapply(mcse, function(x) min(which(x < mc_error_target))))
    # Select the sufficient sample size, with an upper limit of 100
    repetitions <- min(c(nsim[n_ind], 100), na.rm = TRUE)

    cat(
      "The simulation run will be perfomed with", repetitions,
      "repetitions per cell.\n\n"
    )
  }
}


library(pomp)
library(circumstance)
library(doFuture)
library(tidyverse)
plan(multicore)

exp_growth_ode <- Csnippet("
  double dW = rnorm(0, sqrt(dt));
  Y += r * (a - Y) * dt + sigma_dyn * dW;
")

rmeas <- Csnippet("
  Y_obs = rnorm(Y, sigma_meas);
")

init <- Csnippet("
  Y = -2;
")

exp_growth_pomp <- pomp(
  data = data.frame(time = 1:200, Y_obs = NA),
  times = "time",
  t0 = 0,
  rprocess = euler(exp_growth_ode, delta.t = (50 / 7) / 234),
  rinit = init,
  rmeasure = rmeas,
  statenames = "Y",
  paramnames = c("r", "a", "sigma_dyn", "sigma_meas")
)

df <- simulate(
  exp_growth_pomp,
  params = c(r = 0.02, a = 2, sigma_dyn = sqrt(1), sigma_meas = sqrt(1)),
  format = "data.frame"
)

plot(df$time, df$Y_obs)
lines(df$time, df$Y)

exp_growth_pomp@data <- array(df$Y_obs,
  dim = c(1, 200),
  dimnames = list("Y_obs")
)

dmeas <- Csnippet("
  lik = dnorm(Y_obs, Y, sigma_meas, give_log);
")

pf_exp <- pfilter(
  exp_growth_pomp,
  Np = 10000,
  dmeasure = dmeas,
  params = c(r = 0.02, a = 2, sigma_dyn = sqrt(1), sigma_meas = sqrt(1)),
  paramnames = c("sigma_meas"),
  statenames = c("Y")
)

logLik(pf_exp)
plot(pf_exp)
as(pf_exp, "data.frame")

guesses <- sobol_design(
  lower = c(r = 0, a = 0, sigma_dyn = 0, sigma_meas = 0),
  upper = c(r = 5, a = 100, sigma_dyn = 100, sigma_meas = 100),
  nseq = 5
)

plot(guesses, pch = 16)

m_exp <- mif2(
  exp_growth_pomp,
  params = guesses[1, ],
  Np = 1000,
  Nmif = 20,
  dmeasure = dmeas,
  partrans = parameter_trans(log = c("r", "a", "sigma_dyn", "sigma_meas")),
  rw.sd = rw_sd(r = 0.02, a = 0.02, sigma_dyn = 0.02, sigma_meas = 0.02),
  cooling.fraction.50 = 0.5,
  paramnames = c("r", "a", "sigma_dyn", "sigma_meas"),
  statenames = c("Y")
)
plot(m_exp)

logmeanexp(
  logLik(circumstance::pfilter(m_exp, Nrep = 5)),
  se = TRUE, ess = TRUE
)

mifs <- circumstance::mif2(
  m_exp,
  starts = guesses,
  Nmif = 30,
  Np = 2000
)
plot(mifs)

mifs |>
  circumstance::pfilter(Nrep = 5) |>
  logLik() |>
  melt() |>
  separate(name, into = c(".id", "rep")) |>
  group_by(.id) |>
  reframe(melt(logmeanexp(value, se = TRUE))) |>
  ungroup() |>
  bind_rows(
    mifs |> coef() |> melt()
  ) |>
  pivot_wider() |>
  rename(
    loglik = est,
    loglik.se = se
  ) -> estimates
estimates


estimates |>
  bind_rows(guesses) |>
  filter(is.na(loglik) | loglik>max(loglik,na.rm=TRUE)-30) |>
  mutate(col=if_else(is.na(loglik),"#99999955","#ff0000ff")) |>
  {
    \(dat) pairs(~loglik+r+sigma+K+N_0,data=dat,col=dat$col,pch=16)
  }()

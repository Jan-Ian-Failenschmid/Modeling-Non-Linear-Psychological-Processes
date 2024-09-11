#' ----------------------------------------------------------------------------#
#' Title: Simulation Main File                                                 #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 10-04-2024                                                    #
#' -----                                                                       #
#' Last Modified: 11-09-2024                                                   #
#' Last Modified: 11-09-2024                                                   #
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
library(mgcv)
library(cmdstanr)
library(dynr)
library(nprobust)
library(data.table)
library(ggdist)
library(future.apply)

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
papaja::r_refs("./R_bib.bib", append = FALSE)

### Simulation parameters ------------------------------------------------------
## Generative Models
exp_growth <- new("gen_model",
  model_name = "exp_growth",
  model_type = "DE",
  # Set dyn_er to 0 and overwrite during the simulation
  time = 200,
  pars = list(yr = 0.02, ya = 2, dyn_er = 1),
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
  time = 200,
  pars = list(k = 4.3, r = 0.04, dyn_er = 1),
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
  time = 200,
  # Set dyn_er to 0 and overwrite during the simulation
  delta = (50 / 7) / 234, # Integer divides 3, 6, and 9 obs. per day
  stepsize = (50 / 7) / 3, # Gets overwritten anyway
  pars = list(a = -5, dyn_er = 1, omega = -(2 * pi / 50)^2),
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
  time = 200,
  # Set dyn_er to 0 and overwrite during the simulation
  delta = (50 / 7) / 234, # Integer divides 3, 6, and 9 obs. per day
  stepsize = (50 / 7) / 3, # Gets overwritten anyway
  pars = list(k = 0.01, c = 0.1, dyn_er = 1),
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

simple <- new("method_simple",
  method_name = "simple"
)

poly <- new("method_poly",
  method_name = "poly"
)

### Run simulation -------------------------------------------------------------
repetitions <- 30 # Number of repetitions in the pilot sample
mc_error_target <- 0.1 # Desired monte carlo error
for (run in c("pilot")) {
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
      # method_list = list(locpol, gp, gam, dynm, simple, poly),
      method_list = list(gam, dynm),
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
    rm(sim) # Remove sim from working environment
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
    repetitions <- min(c(nsim[n_ind], 5), na.rm = TRUE)

    cat(
      "The simulation run will be perfomed with", repetitions,
      "repetitions per cell.\n\n"
    )
  }
}

#' ----------------------------------------------------------------------------#
#' Title: Simulation Main File                                                 #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 10-04-2024                                                    #
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

poly_orth <- new("method_poly_orth",
  method_name = "poly_orth"
)

### Run simulation -------------------------------------------------------------
for (run in c("simulation")) {
  # Set seed
  if (run == "pilot") {
    repetitions <- 30
    set.seed(12345)
  } else if (run == "simulation") {
    repetitions <- 100
    set.seed(54321)
  }

  # Run simulation
  cat("Running", run, "---------------------------------------------------\n\n")
  system.time({
    sim <- simulate(
      gen_model_list = list(
        exp_growth, log_growth, damped_oscillator, cusp_catastrophe
      ),
      method_list = list(locpol, gp, gam, dynm, simple, poly, poly_orth),
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
}

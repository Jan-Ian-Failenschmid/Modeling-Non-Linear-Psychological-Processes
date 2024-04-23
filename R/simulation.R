#' ----------------------------------------------------------------------------#
#' Title: Simulation Main File                                                 #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 10-04-2024                                                    #
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
  dynr, # State-space modelling
  nprobust, # Local polynomial estimator
  parallel # Parallel computing
)

# Load functions
invisible(sapply(
  c(paste0("./R/helper/", dir(path = "./R/helper"))),
  source
))

### Simulation parameters ------------------------------------------------------
## Generative Models
latent_change <- new("gen_model",
  model_name = "latent_change",
  model_type = "DE",
  pars = list(yr = 0.02, ya = 2, dyn_er = 0),
  delta = 1 / 32,
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

cusp_catastrophe <- new("gen_model",
  model_name = "cusp_catastrophe",
  model_type = "DE",
  delta = 1 / 32,
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
  delta = 1 / 32,
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
system.time({
  sim <- simulate(
    gen_model_list = list(
      latent_change, log_growth, damped_oscillator, cusp_catastrophe
    ),
    method_list <- list(locpol, gp, gam, ssm),
    conditions = list(
      time = c(50, 100),
      stepsize = c(0.5, 1),
      dyn_er = c(0.1, 0.25)
    ),
    repetitions = 1
  )
})

# Check object size of sim grid
cat(utils::object.size(sim)[1] / (1e6), "mb")

### Save simulated data --------------------------------------------------------
sim_time <- format(Sys.time(), "%d_%m_%Y_%H_%M")
save(sim, file = paste0("./R/data/simulated_data_", sim_time, ".Rdata"))

### Extract and save results ---------------------------------------------------
res <- extract_results(sim)
save(res, file = paste0("./R/data/simulation_results_", sim_time, ".Rdata"))

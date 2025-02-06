#' ----------------------------------------------------------------------------#
#' Title:                                                                      #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 03-02-2025                                                    #
#' -----                                                                       #
#' Last Modified: 06-02-2025                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2025 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#

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

# Restart simulation
cat("Restarting ----------------------------------------------------------\n\n")
system.time({
  sim <- restart_simulation(out_dir = out_dir)
})

load(paste0(out_dir, "/temp.Rdata"))
sim <- sim_grid
# Save results
# Check object size of sim grid
cat(utils::object.size(sim)[1] / (1e6), "mb")

# Save simulation
sim_time <- format(Sys.time(), "%d_%m_%Y_%H_%M")
save(sim, file = paste0(out_dir, "/simulation_data_", sim_time, ".Rdata"))

# Extract and save results
res <- extract_results(sim)
save(res, file = paste0(out_dir, "/simulation_results_", sim_time, ".Rdata"))

# Combined data sets
load(paste0(out_dir, "/simulation_data_27_09_2024_00_32.Rdata"))
old_sim <- sim

load(paste0(out_dir, "/simulation_data_06_02_2025_19_31.Rdata"))
sim$method[[1]][[5]]

for (i in seq_len(nrow(old_sim))) {
  old_sim$method[[i]][[2]] <- sim$method[[i]][[2]]
  print(i)
}

mapply(function(x, y) {
  identical(x[[2]], y[[2]])
}, x = old_sim$method, y = sim$method)

old_sim$method[[4801]][[4]]
sim$method[[4801]][[4]]

# Save combined data set
sim <- old_sim
sim_time <- format(Sys.time(), "%d_%m_%Y_%H_%M")
save(sim, file = paste0(out_dir, "/combined_data_", sim_time, ".Rdata"))

# Extract and save results
res <- extract_results(sim)
save(res, file = paste0(out_dir, "/combined_results_", sim_time, ".Rdata"))

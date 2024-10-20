#' ----------------------------------------------------------------------------#
#' Title:                                                                      #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 20-10-2024                                                    #
#' -----                                                                       #
#' Last Modified: 20-10-2024                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#

server <- function(input, output) {
  # Read in simulation data and results
  sim <- load_sim_data()
  load("./www/data/simulation_results_27_09_2024_00_32.Rdata")
  res <- as.data.table(res)
  res <- res[!method %in% c("simple", "poly_orth"), ]
}

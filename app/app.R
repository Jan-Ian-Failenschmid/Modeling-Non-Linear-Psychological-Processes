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

# Load packages
library(shiny)

# Source helper files
source("./www/classes.R")
source("./www/functions.R")

# Source UI
source("ui.R")

# Source server
source("server.R")

# Create app
shinyApp(ui = ui, server = server)

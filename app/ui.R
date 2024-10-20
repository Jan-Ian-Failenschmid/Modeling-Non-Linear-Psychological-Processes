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

ui <- fluidPage(
  # Header
  headerPanel(h1("Modeling Non-Linear Psychological Processes")),
  # Tabset to switch in between pages
  tabsetPanel(
    id = "switcher",
    # First panel contains the introduction
    tabPanel(
      "Introduction",
      h4("This application contains the online supplementary material to the
      article .")
    ),
    # Second panel contains the smoothing plots
    tabPanel(
      "Smoothing Plots",
      sidebarLayout(
        sidebarPanel(
          fluidRow(
            checkboxGroupInput("smooth_methods",
              "Which methods would you like to include?",
              choices = list(
                "Local Polynomial Regression",
                "Generalized Additive Model",
                "Gaussian Process",
                "Polynomial Regression",
                "Linear Regression",
                "Parametric Models"
              ), selected = list(
                "Local Polynomial Regression",
                "Generalized Additive Model",
                "Gaussian Process",
                "Polynomial Regression",
                "Linear Regression",
                "Parametric Models"
              )
            )
          ),
          fluidRow(
            checkboxGroupInput("smooth_processes",
              "Which processes would you like to include?",
              choices = list(
                "Exponential Growth",
                "Logistic Growth",
                "Cusp Catastrophe",
                "Damped Oscillator"
              ), selected = list(
                "Exponential Growth",
                "Logistic Growth",
                "Cusp Catastrophe",
                "Damped Oscillator"
              )
            )
          ),
          fluidRow(
            numericInput("smooth_numbers",
              "Which data set would you like to plot?",
              value = 1,
              min = 1,
              max = 30
            )
          )
        ),
        mainPanel()
      )
    ),
    tabPanel(
      "Results",
      sidebarLayout(
        sidebarPanel(
          fluidRow(
            selectInput("results_dependent",
              "Which outcome variable would you like to look at?",
              choices = list("MSE", "GCV", "CI-Coverage"),
              selected = "MSE"
            )
          ),
          fluidRow(
            checkboxGroupInput("results_independent",
              "Which predictors would you like to use?",
              choices = list(
                "Method",
                "Process",
                "Dynamic Error Variance",
                "Measurment Period",
                "Measurment Frequency"
              ), selected = list(
                "Method",
                "Process"
              )
            )
          )
        ),
        mainPanel()
      )
    )
  )
)

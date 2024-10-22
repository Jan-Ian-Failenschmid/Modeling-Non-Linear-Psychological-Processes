#' ----------------------------------------------------------------------------#
#' Title:                                                                      #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 20-10-2024                                                    #
#' -----                                                                       #
#' Last Modified: 22-10-2024                                                   #
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
library(data.table)
library(ggplot2)
library(ggdist)
library(papaja)

# Source helper files
invisible(sapply(
  c(paste0("./helper/", dir(path = "./helper"))),
  source
))
# Define necessary methods
for (class_name in c(
  "method_poly", "method_gp",
  "method_poly_orth", "method_simple", "method_dynm"
)) {
  setClass(
    class_name,
    contains = "method"
  )
}

# Read-in simulation data
load("data/simulation_data_27_09_2024_00_32.Rdata")

# Read-in simulation results
load("./data/simulation_results_27_09_2024_00_32.Rdata")
res <- as.data.table(res)


res[, SP := ifelse(time == 100, 0.5, 1)]
res[, DEV := dyn_er^2]
res[, SF := (1 / (stepsize * (7 / 50))) / 3]


res$model[res$model == "exp_growth"] <- "Exponential Growth"
res$model[res$model == "log_growth"] <- "Logistic Growth"
res$model[res$model == "damped_oscillator"] <- "Damped Oscillator"
res$model[res$model == "cusp_catastrophe"] <- "Cusp Catastrophe"

res$method[res$method == "locpol"] <- "Local Polynomial Regression"
res$method[res$method == "gp"] <- "Gaussian Process Regression"
res$method[res$method == "gam"] <- "Generalized Additive Modeling"
res$method[res$method == "dynm"] <- "Parametric Modeling"
res$method[res$method == "simple"] <- "Linear Regression"
res$method[res$method == "poly"] <- "Global Polynomial Regression"
res$method[res$method == "poly_orth"] <- "Orthogonal Polynomial Regression"

res$model <- factor(res$model,
  levels = c(
    "Exponential Growth", "Logistic Growth",
    "Damped Oscillator", "Cusp Catastrophe"
  )
)
contrasts(res$model) <- contr.sum(levels(res$model))

res$method <- factor(res$method,
  levels = c(
    "Local Polynomial Regression", "Gaussian Process Regression",
    "Generalized Additive Modeling",
    "Global Polynomial Regression",
    "Orthogonal Polynomial Regression",
    "Linear Regression", "Parametric Modeling"
  )
)
contrasts(res$method) <- contr.sum(levels(res$method))

res$SP <- factor(res$SP)
contrasts(res$SP) <- contr.sum(levels(res$SP))

res$SF <- factor(res$SF)
contrasts(res$SF) <- contr.sum(levels(res$SF))

res$DEV <- factor(res$DEV,
  levels = c(.5, 1, 2),
  labels = c(".5", "1", "2")
)
contrasts(res$DEV) <- contr.sum(levels(res$DEV))

# Method and process lists
method_list <- c(
  "Local Polynomial Regression",
  "Gaussian Process Regression",
  "Generalized Additive Modeling",
  "Parametric Modeling",
  "Linear Regression",
  "Global Polynomial Regression",
  "Orthogonal Polynomial Regression"
)

process_list <- c(
  "Exponential Growth",
  "Logistic Growth",
  "Damped Oscillator",
  "Cusp Catastrophe"
)

# UI ---------------------------------------------------------------------------
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
              choiceNames = list(
                "Local Polynomial Regression",
                "Gaussian Process Regression",
                "Generalized Additive Modeling",
                "Global Polynomial Regression",
                "Orthogonal Polynomial Regression",
                "Linear Regression",
                "Parametric Modeling"
              ),
              choiceValues = list(1, 2, 3, 6, 7, 5, 4),
              selected = list(1, 2, 3, 4, 6)
            )
          ),
          fluidRow(
            checkboxGroupInput("smooth_processes",
              "Which processes would you like to include?",
              choiceNames = list(
                "Exponential Growth",
                "Logistic Growth",
                "Damped Oscillator",
                "Cusp Catastrophe"
              ),
              choiceValues = list(1, 2, 3, 4),
              selected = list(1, 2, 3, 4)
            )
          ),
          fluidRow(
            h4("Which condition would you like to plot?"),
            selectInput(
              "smooth_dyn_var",
              "Dynamic Error Variance:",
              choices = list(
                0.5, 1, 2
              ), selected = 0.5
            ),
            selectInput(
              "smooth_period",
              "Sampling Period:",
              choices = list(
                "Full", "Half"
              ), selected = "Full"
            ),
            selectInput(
              "smooth_frequency",
              "Sampling Frequency:",
              choices = list(
                1, 2, 3
              ), selected = 3
            )
          ),
          fluidRow(
            numericInput("smooth_number",
              "Which data set would you like to plot?",
              value = 1,
              min = 1,
              max = 100
            )
          ),
          fluidRow(
            h4("Adjust graph size:"),
            numericInput("smooth_graph_hight",
              "Graph hight:",
              value = 600,
              min = 400,
              max = 2000
            ),
            numericInput("smooth_graph_width",
              "Graph width:",
              value = 800,
              min = 400,
              max = 2000
            )
          ),
          width = 2
        ),
        mainPanel(
          fluidRow(
            h1("Smoothing Plots:")
          ),
          fluidRow(
            plotOutput("smoothing_plot", width = "auto", height = "auto")
          ),
        )
      )
    ),
    tabPanel(
      "Results",
      sidebarLayout(
        sidebarPanel(
          fluidRow(
            selectInput("results_style",
              "Which graph styile would you like to have?",
              choices = list(
                "Complete" = "all",
                "Average" = "mean",
                "Missing" = "missing"
              ),
              selected = "mean"
            )
          ),
          fluidRow(
            selectInput("results_dependent",
              "Which outcome variable would you like to look at?",
              choices = list(
                "MSE" = "mse",
                "GCV" = "gcv",
                "CI-Coverage" = "ci_coverage"
              ),
              selected = "MSE"
            )
          ),
          fluidRow(
            checkboxGroupInput("results_independent",
              "Which predictors would you like to use?",
              choices = list(
                "Dynamic Error Variance",
                "Sampling Period",
                "Sampling Frequency"
              )
            )
          ),
          fluidRow(
            h4("Select which levels of each factor you would like to include:"),
            checkboxGroupInput("results_methods",
              "Methods:",
              choices = list(
                "Local Polynomial Regression",
                "Gaussian Process Regression",
                "Generalized Additive Modeling",
                "Global Polynomial Regression",
                "Orthogonal Polynomial Regression",
                "Linear Regression",
                "Parametric Modeling"
              ),
              selected = list(
                "Local Polynomial Regression",
                "Gaussian Process Regression",
                "Generalized Additive Modeling",
                "Global Polynomial Regression",
                "Parametric Modeling"
              )
            ),
            checkboxGroupInput("results_processes",
              "Processes:",
              choices = list(
                "Exponential Growth",
                "Logistic Growth",
                "Damped Oscillator",
                "Cusp Catastrophe"
              ),
              selected = list(
                "Exponential Growth",
                "Logistic Growth",
                "Damped Oscillator",
                "Cusp Catastrophe"
              )
            ),
            uiOutput(
              "dynamic_error_levles"
            ),
            uiOutput(
              "sampling_period_levles"
            ),
            uiOutput(
              "sampling_frequency_levles"
            )
          ),
          width = 2
        ),
        mainPanel(
          plotOutput("results_plot", width = "auto", height = "auto")
        )
      )
    ),
    tabPanel(
      "ANOVA",
      sidebarLayout(
        sidebarPanel(
          fluidRow(
            selectInput("anova_dependent",
              "Which outcome variable would you like to look at?",
              choices = list(
                "MSE" = "mse",
                "GCV" = "gcv"
              ),
              selected = "MSE"
            )
          ),
          fluidRow(
            checkboxGroupInput("anova_independent",
              "Which predictors would you like to use?",
              choices = list(
                "Dynamic Error Variance",
                "Sampling Period",
                "Sampling Frequency"
              )
            )
          ),
          fluidRow(
            h4("Select which levels of each factor you would like to include:"),
            checkboxGroupInput("anova_methods",
              "Methods:",
              choices = list(
                "Local Polynomial Regression",
                "Gaussian Process Regression",
                "Generalized Additive Modeling",
                "Global Polynomial Regression",
                "Orthogonal Polynomial Regression",
                "Linear Regression"
              ),
              selected = list(
                "Local Polynomial Regression",
                "Gaussian Process Regression",
                "Generalized Additive Modeling",
                "Global Polynomial Regression",
                "Parametric Modeling"
              )
            ),
            checkboxGroupInput("anova_processes",
              "Processes:",
              choices = list(
                "Exponential Growth",
                "Logistic Growth",
                "Damped Oscillator",
                "Cusp Catastrophe"
              ),
              selected = list(
                "Exponential Growth",
                "Logistic Growth",
                "Damped Oscillator",
                "Cusp Catastrophe"
              )
            ),
            uiOutput(
              "anova_dynamic_error_levles"
            ),
            uiOutput(
              "anova_sampling_period_levles"
            ),
            uiOutput(
              "anova_sampling_frequency_levles"
            )
          ),
          width = 2
        ),
        mainPanel(
          verbatimTextOutput("test")
        )
      )
    ),
    tabPanel(
      "Data Generation and Model Fitting",
      sidebarLayout(
        sidebarPanel(
          tabsetPanel(
            id = "gen_and_switch",
            tabPanel(
              "Data Generation",
              fluidRow(
                h4("Seed:"),
                checkboxInput("dg_set_seed", "Set Seed"),
                uiOutput(
                  "set_seed_input"
                )
              ),
              fluidRow(
                h4("Process:"),
                selectInput("dg_process_choice",
                  "From which process would you like to generate data?",
                  choices = list(
                    "Exponential Growth",
                    "Logistic Growth",
                    "Damped Oscillator",
                    "Cusp Catastrophe"
                  ),
                  selected = "Exponential Growth"
                )
              ),
              fluidRow(
                h4("Model Parameters:")
              ),
              uiOutput(
                "data_generation_input"
              ),
              fluidRow(
                numericInput("sigma_eta", "Dynamic Error Variance
                \\(\\sigma_\\eta\\):",
                  value = 1
                ),
                numericInput("sigma_epsilon", withMathJax(
                  "Measurment Error Variance \\(\\sigma_\\epsilon\\):"
                ),
                value = 1
                )
              ),
              fluidRow(
                h4("Simulation Parameters:"),
                numericInput("t_max", "Time Max \\(t_{max}\\):",
                  value = 200, min = 0
                ),
                numericInput("t_step", "Time Step Lenght \\(\\Delta t\\):",
                  value = 1, min = 0
                ),
              ),
              fluidRow(
                h4("Euler-Maruyama Method"),
                numericInput("subdiv", "Subdivisions per time steps
                \\(\\delta\\):",
                  value = 200, min = 1
                ),
                h4("Simulate Data"),
                actionButton("sim_data", "Simulate Data"),
                h5("This may take a while.")
              )
            ),
            tabPanel(
              "Model Fitting",
              fluidRow(
                h4("Which methods would you like to fit?"),
                h4("Local Polynomial Regression:"),
                checkboxInput("fit_lpr", "Fit Local Polynomial Regression"),
                uiOutput("locpol_settings_ui"),
                h4("Generalized Additive Model:"),
                checkboxInput("fit_gam", "Fit Generalized Additive Model"),
                uiOutput("gam_settings_ui"),
                h4("Fit models:"),
                actionButton("fit_models", "Fit Models"),
                h5("This may take a while.")
              )
            )
          ),
          width = 2
        ),
        mainPanel(
          fluidRow(
            h2("Data Generation:"),
            column(
              width = 7,
              h3("Model Equations:"),
              uiOutput(
                "model_equations"
              ),
              h3("Euler-Maruyama Approximation:"),
              uiOutput(
                "em_equations"
              ),
              h3("Measurement Equations:"),
              uiOutput(
                "meas_equations"
              ),
              offset = 0
            ),
            column(
              width = 5,
              uiOutput("tsm_model_header"),
              plotOutput("tsm_model", width = "auto", height = "auto"),
              offset = 0
            )
          ),
          uiOutput("fit_ui_header"),
          uiOutput("gam_ui_pop"),
          uiOutput("locpol_ui_pop")
        )
      )
    )
  )
)

# Server -----------------------------------------------------------------------
server <- function(input, output) {
  ## Smoothing Plots -----------------------------------------------------------
  # Find data set indicators
  ilustr <- reactive(sim[, model_name := sapply(gen_model, function(x) {
    x@model_name
  })][,
    .I[time == ifelse(input$smooth_period == "Half", 100, 200) &
      stepsize == unique(stepsize)[as.numeric(input$smooth_frequency)] &
      dyn_er == sqrt(as.numeric(input$smooth_dyn_var))],
    by = model_name
  ][, .SD[input$smooth_number], by = model_name])

  # Render Plot
  output$smoothing_plot <- renderPlot(
    {
      par(mfrow = c(
        length(input$smooth_processes),
        length(input$smooth_methods)
      ))

      for (j in as.numeric(input$smooth_processes)) {
        for (i in as.numeric(input$smooth_methods)) {
          method <- sim$method[[ilustr()$V1[j]]][[i]]

          converged <- method@converged
          if (converged) {
            plot(method,
              sim = sim,
              axes = FALSE, xlab = " ", ylab = ""
            )
          } else {
            plot_dat <- sim$dat[[ilustr()$V1[j]]]
            plot(
              x = plot_dat$time, y = plot_dat$y_obs,
              axes = FALSE, xlab = "", ylab = ""
            )
            lines(x = plot_dat$time, y = plot_dat$y)
          }

          axis(2, at = seq(-10, 10, 2), cex.axis = 2.5)

          axis(1,
            at = c(0, 50, 100, 150, 200),
            labels = c("", "First Half", "", "Sec. Half", ""), cex.axis = 2.5,
            padj = 1
          )

          if (j == 1) {
            title(method_list[i])
          } else if (j == length(input$smooth_processes)) {
            title(xlab = "Time", main = NULL, line = 5)
          }
          if (i == min(as.numeric(input$smooth_methods))) {
            title(ylab = process_list[j], main = NULL)
          } else {
            title(main = NULL)
          }
        }
      }
    },
    width = function() {
      input$smooth_graph_width * length(input$smooth_methods)
    },
    height = function() {
      input$smooth_graph_hight * length(input$smooth_processes)
    }
  )

  ## Results plot --------------------------------------------------------------
  # Select input levels
  output$dynamic_error_levles <- renderUI({
    if ("Dynamic Error Variance" %in% input$results_independent) {
      checkboxGroupInput("results_dyn_var",
        "Dynamic Error Variance:",
        choices = list(
          ".5",
          "1",
          "2"
        ),
        selected = list(
          ".5",
          "1",
          "2"
        )
      )
    }
  })

  output$sampling_period_levles <- renderUI({
    if ("Sampling Period" %in% input$results_independent) {
      checkboxGroupInput("results_sampling_period",
        "Sampling Period:",
        choiceNames = list(
          "Half",
          "Full"
        ),
        choiceValues = list(
          0.5, 1
        ),
        selected = list(
          0.5, 1
        )
      )
    }
  })

  output$sampling_frequency_levles <- renderUI({
    if ("Sampling Frequency" %in% input$results_independent) {
      checkboxGroupInput("results_sampling_frequency",
        "Sampling Frequency:",
        choices = list(
          "1",
          "2",
          "3"
        ),
        selected = list(
          "1",
          "2",
          "3"
        )
      )
    }
  })

  # Render plot
  output$results_plot <- renderPlot(
    {
      res_plot <- res[method %in% input$results_methods &
        model %in% input$results_processes, ]
      if ("Dynamic Error Variance" %in% input$results_independent) {
        DEV <- "DEV"
        res_plot <- res_plot[DEV %in% input$results_dyn_var, ]
      } else {
        DEV <- NULL
      }
      if ("Sampling Period" %in% input$results_independent) {
        SP <- "SP"
        res_plot <- res_plot[SP %in% input$results_sampling_period, ]
      } else {
        SP <- NULL
      }
      if ("Sampling Frequency" %in% input$results_independent) {
        SF <- "SF"
        res_plot <- res_plot[SF %in% input$results_sampling_frequency, ]
      } else {
        SF <- NULL
      }

      p <- plot_results(
        res = res_plot,
        input$results_dependent,
        input$results_style, SP, SF, DEV
      )

      if (input$results_dependent == "mse") {
        y_lab_str <- "Mean MSE"
      } else if (input$results_dependent == "gcv") {
        y_lab_str <- "Mean GCV"
      } else if (input$results_dependent == "ci_coverage") {
        y_lab_str <- "Mean CI-Coverage"
        p <- p + annotate("rect",
          xmin = -Inf, xmax = Inf, ymin = 0.89, ymax = 1,
          alpha = 0.3
        )
      }

      p <- p +
        theme_apa() +
        theme(
          text = element_text(size = 26)
        ) + labs(
          title = "",
          y = y_lab_str, fill = "Process", color = "Process",
          x = "Simulation Conditions"
        )

      p
    },
    width = 1800,
    height = 1000
  )

  ## ANOVA ---------------------------------------------------------------------
  output$anova_dynamic_error_levles <- renderUI({
    if ("Dynamic Error Variance" %in% input$anova_independent) {
      checkboxGroupInput("results_dyn_var",
        "Dynamic Error Variance:",
        choices = list(
          ".5",
          "1",
          "2"
        ),
        selected = list(
          ".5",
          "1",
          "2"
        )
      )
    }
  })

  output$anova_sampling_period_levles <- renderUI({
    if ("Sampling Period" %in% input$anova_independent) {
      checkboxGroupInput("results_sampling_period",
        "Sampling Period:",
        choiceNames = list(
          "Half",
          "Full"
        ),
        choiceValues = list(
          0.5, 1
        ),
        selected = list(
          0.5, 1
        )
      )
    }
  })

  output$anova_sampling_frequency_levles <- renderUI({
    if ("Sampling Frequency" %in% input$anova_independent) {
      checkboxGroupInput("results_sampling_frequency",
        "Sampling Frequency:",
        choices = list(
          "1",
          "2",
          "3"
        ),
        selected = list(
          "1",
          "2",
          "3"
        )
      )
    }
  })


  output$test <- renderPrint({
    contrasts(res$method)
  })

  ## Data Generation and model fitting -----------------------------------------
  ### Model equations ----------------------------------------------------------
  output$set_seed_input <- renderUI({
    if (input$dg_set_seed) {
      numericInput("dg_seed", "Seed Value:", value = 1, min = 1)
    }
  })

  output$data_generation_input <- renderUI({
    if (input$dg_process_choice == "Exponential Growth") {
      fluidRow(
        numericInput("y_nod", withMathJax("Starting Value \\(y_0\\):"),
          value = -2
        ),
        numericInput("r", withMathJax("Growth Rate \\(\\theta\\):"),
          value = 0.02, min = 0
        ),
        numericInput("a", "Asymptote \\(\\mu\\):", value = 2),
      )
    } else if (input$dg_process_choice == "Logistic Growth") {
      fluidRow(
        numericInput("y_nod", withMathJax("Starting Value \\(y_0\\):"),
          value = 0.01, min = 0
        ),
        numericInput("r", withMathJax("Growth Rate \\(\\theta\\):"),
          value = 0.04, min = 0
        ),
        numericInput("a", "Asymptote \\(\\mu\\):", value = 4.3, min = 0),
      )
    } else if (input$dg_process_choice == "Damped Oscillator") {
      fluidRow(
        numericInput("y_nod", withMathJax("Starting Value \\(y_0\\):"),
          value = 2
        ),
        numericInput("v_nod", withMathJax("Starting Value \\(v_0\\):"),
          value = 0
        ),
        numericInput("c", withMathJax("Growth Rate \\(c\\):"),
          value = 0.1, min = 0
        ),
        numericInput("k", "Asymptote \\(k\\):",
          value = 0.01, min = 0
        ),
      )
    } else if (input$dg_process_choice == "Cusp Catastrophe") {
      fluidRow(
        numericInput("y_nod", withMathJax("Starting Value \\(y_0\\):"),
          value = 1.9
        ),
        numericInput("b_nod", withMathJax("Starting Value \\(b_0\\):"),
          value = -10
        ),
        numericInput("v_nod", withMathJax("Starting Value \\(v_0\\):"),
          value = 0
        ),
        numericInput("a", withMathJax("Growth Rate \\(a\\):"),
          value = -5, min = 0
        ),
        numericInput("omega", "Asymptote \\(\\omega\\):",
          value = -0.016, min = 0
        ),
      )
    }
  })

  output$model_equations <- renderUI({
    if (input$dg_process_choice == "Exponential Growth") {
      withMathJax(h4(paste0(
        "$$
      \\begin{aligned}
        &dy = \\theta  (\\mu - y) \\, dt + \\sigma_\\eta \\, dW_t \\\\
        \\\\
        &dy = ",
        input$r, " (",
        input$a, " - y) \\, dt + ",
        input$sigma_eta, " \\, dW_t
      \\end{aligned}
      $$"
      )))
    } else if (input$dg_process_choice == "Logistic Growth") {
      withMathJax(h4(paste0(
        "$$
      \\begin{aligned}
        &dy = \\theta y (1 - \\frac{y}{\\mu}) \\, dt +
        \\sigma_\\eta \\, dW_t \\\\
        \\\\
        &dy = ",
        input$r, " y (1 - \\frac{y}{",
        input$a, "}) \\, dt + ",
        input$sigma_eta, " \\, dW_t
      \\end{aligned}
      $$"
      )))
    } else if (input$dg_process_choice == "Damped Oscillator") {
      withMathJax(h4(paste0(
        "$$
      \\begin{aligned}
        &dy = v \\, dt \\\\
        &dv = (-2 k v - c^2 y) \\, dt + \\sigma_\\eta \\, dW_t \\\\
        \\\\
        &dy = v \\, dt \\\\
        &dv = (-2 * ",
        input$k, " v - ",
        input$c, "^2 y) \\, dt + ",
        input$sigma_eta, " \\, dW_t
      \\end{aligned}
      $$"
      )))
    } else if (input$dg_process_choice == "Cusp Catastrophe") {
      withMathJax(h4(paste0(
        "$$
      \\begin{aligned}
        &dy = -(4  y^3 + 2  a  y + b) \\, dt + \\sigma_\\eta \\, dW_t\\\\
        &db = v \\, dt \\\\
        &dv = \\omega b \\, dt \\\\
        \\\\
        &dy = -(4 y^3 + 2 * ",
        input$a, " * y + b) \\, dt + ",
        input$sigma_eta, " \\, dW_t\\\\
        &db = v \\, dt \\\\
        &dv = ",
        input$omega, " b \\, dt
      \\end{aligned}
      $$"
      )))
    }
  })

  output$em_equations <- renderUI({
    if (input$dg_process_choice == "Exponential Growth") {
      withMathJax(h4(paste0(
        "$$
      \\begin{aligned}
        &y_{t+\\frac{\\Delta t}{\\delta}} = y_t +
        \\frac{\\Delta t}{\\delta} (\\theta  (\\mu - y_t)) +
        \\sigma_\\eta \\eta_t \\\\
        &\\eta_t \\sim N(0, \\frac{\\Delta t}{\\delta}) \\\\
        \\\\

        &y_{t+ \\frac{", input$t_step, "}{", input$subdiv, "}} = y_t +
        \\frac{", input$t_step, "}{", input$subdiv, "} (",
        input$r, "  (", input$a, " - y_t)) + ",
        input$sigma_eta, " \\eta_t \\\\
        &\\eta_t \\sim N(0,
        \\frac{", input$t_step, "}{", input$subdiv, "})
      \\end{aligned}
      $$"
      )))
    } else if (input$dg_process_choice == "Logistic Growth") {
      withMathJax(h4(paste0(
        "$$
      \\begin{aligned}
        &y_{t+\\frac{\\Delta t}{\\delta}} = y_t +
        \\frac{\\Delta t}{\\delta} (\\theta y (1 - \\frac{y}{\\mu})) +
        \\sigma_\\eta \\eta_t \\\\
        &\\eta_t \\sim N(0, \\frac{\\Delta t}{\\delta}) \\\\
        \\\\

        &y_{t+ \\frac{", input$t_step, "}{", input$subdiv, "}} = y_t +
        \\frac{", input$t_step, "}{", input$subdiv, "}(",
        input$r, " y (1 - \\frac{y_t}{", input$a, "})) + ",
        input$sigma_eta, " \\eta_t \\\\
        &\\eta_t \\sim N(0, \\frac{", input$t_step, "}{", input$subdiv, "})
      \\end{aligned}
      $$"
      )))
    } else if (input$dg_process_choice == "Damped Oscillator") {
      withMathJax(h4(paste0(
        "$$
      \\begin{aligned}
        &y_{t + \\frac{\\Delta t}{\\delta}} = y_t +
        \\frac{\\Delta t}{\\delta} v_t +
        \\sigma_\\eta \\eta_t \\\\\
        &v_{t+\\frac{\\Delta t}{\\delta}} = v_t +
        \\frac{\\Delta t}{\\delta} (-2 k v_t - c^2 y_t) \\\\
        &\\eta_t \\sim N(0, \\frac{\\Delta t}{\\delta}) \\\\
        \\\\
        &y_{t+\\frac{", input$t_step, "}{", input$subdiv, "}} = y_t +
        \\frac{", input$t_step, "}{", input$subdiv, "} v_t + ",
        input$sigma_eta, " \\eta_t \\\\
        &v_{t+\\frac{", input$t_step, "}{", input$subdiv, "}} = v_t +
        \\frac{", input$t_step, "}{", input$subdiv, "}(-2 * ",
        input$k, " v_t - ",
        input$c, "^2 y_t) \\\\
        &\\eta_t \\sim N(0, \\frac{", input$t_step, "}{", input$subdiv, "})
      \\end{aligned}
      $$"
      )))
    } else if (input$dg_process_choice == "Cusp Catastrophe") {
      withMathJax(h4(paste0(
        "$$
      \\begin{aligned}
        &h =  \\frac{\\Delta t}{\\delta} \\\\
        &y_{t+\\frac{\\Delta t}{\\delta}} = y_t -
        \\frac{\\Delta t}{\\delta} (4  y_t^3 + 2  a  y + b_t) +
        \\sigma_\\eta \\eta_t \\\\
        &b_{t+\\frac{\\Delta t}{\\delta}} = b_t +
        \\frac{\\Delta t}{\\delta} v_t \\\\
        &v_{t+\\frac{\\Delta t}{\\delta}} = v_t +
        \\frac{\\Delta t}{\\delta} (\\omega b_t) \\\\
        &\\eta_t \\sim N(0, \\frac{\\Delta t}{\\delta}) \\\\
        \\\\
        &",
        input$t_step / input$subdiv,
        " = \\frac{", input$t_step, "}{", input$subdiv, "} \\\\
        &y_{t+\\frac{", input$t_step, "}{", input$subdiv, "}} = y_t -
        \\frac{", input$t_step, "}{", input$subdiv, "}(4 y_t^3 + 2 + ",
        input$a, " + y + b_t) +",
        input$sigma_eta, " \\eta_t \\\\
        &b_{t+\\frac{", input$t_step, "}{", input$subdiv, "}} = b_t +
        \\frac{", input$t_step, "}{", input$subdiv, "}v_t \\\\
        &v_{t+\\frac{", input$t_step, "}{", input$subdiv, "}} = v_t +
        \\frac{", input$t_step, "}{", input$subdiv, "}(",
        input$omega, " b_t) \\\\
        &\\eta_t \\sim N(0, \\frac{", input$t_step, "}{", input$subdiv, "})
      \\end{aligned}
      $$"
      )))
    }
  })

  output$meas_equations <- renderUI({
    withMathJax(h4(paste0(
      "$$
      \\begin{aligned}
        &y^{Obs}_{t} = y_t + \\sigma_\\epsilon \\epsilon_t, \\quad
        \\epsilon_t \\sim N(0, 1) \\\\
        \\\\
        &y^{Obs}_{t} = y_t + ",
      input$sigma_epsilon,
      " \\epsilon_t, \\quad \\epsilon_t \\sim N(0, 1)
      \\end{aligned}
      $$"
    )))
  })

  ### Data generation ----------------------------------------------------------
  gen_tsm <- reactive({
    if (input$dg_process_choice == "Exponential Growth") {
      new("gen_model",
        model_name = "exp_growth",
        model_type = "DE",
        # Set dyn_er to 0 and overwrite during the simulation
        time = input$t_max,
        pars = list(
          yr = input$r,
          ya = input$a,
          dyn_er = sqrt(input$sigma_eta),
          meas_er = sqrt(input$sigma_epsilon),
          y_nod = input$y_nod
        ),
        delta = input$t_step / input$subdiv,
        stepsize = input$t_step, # Gets overwritten anyway
        start = list(
          formula(y ~ y_nod)
        ),
        state_eq = list(
          formula(y ~ yr * ya - yr * y)
        ),
        dynamic_error = list(
          formula(y ~ dyn_er)
        ),
        meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, meas_er)))
      )
    } else if (input$dg_process_choice == "Logistic Growth") {
      new("gen_model",
        model_name = "log_growth",
        model_type = "DE",
        # Set dyn_er to 0 and overwrite during the simulation
        time = input$t_max,
        pars = list(
          k = input$a,
          r = input$r,
          dyn_er = sqrt(input$sigma_eta),
          meas_er = sqrt(input$sigma_epsilon),
          y_nod = input$y_nod
        ),
        delta = input$t_step / input$subdiv,
        stepsize = input$t_step, # Gets overwritten anyway
        start = list(
          formula(y ~ y_nod)
        ),
        state_eq = list(
          formula(y ~ r * y * (1 - (y / k)))
        ),
        meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, meas_er))),
        dynamic_error = list(
          formula(y ~ dyn_er)
        ),
        boundary = function(x) ifelse(x >= 0, x, 0)
      )
    } else if (input$dg_process_choice == "Damped Oscillator") {
      new("gen_model",
        model_name = "damped_oscillator",
        model_type = "DE",
        time = input$t_max,
        # Set dyn_er to 0 and overwrite during the simulation
        delta = input$t_step / input$subdiv,
        stepsize = input$t_step, # Gets overwritten anyway
        pars = list(
          k = input$k,
          c = input$c,
          dyn_er = sqrt(input$sigma_eta),
          meas_er = sqrt(input$sigma_epsilon),
          y_nod = input$y_nod,
          v_nod = input$v_nod
        ),
        start = list(
          formula(y ~ y_nod),
          formula(v ~ v_nod)
        ),
        state_eq = list(
          formula(y ~ v),
          formula(v ~ -2 * k * v - c^2 * y)
        ),
        dynamic_error = list(
          formula(y ~ dyn_er),
          formula(v ~ 0)
        ),
        meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, meas_er)))
      )
    } else if (input$dg_process_choice == "Cusp Catastrophe") {
      new("gen_model",
        model_name = "cusp_catastrophe",
        model_type = "DE",
        time = input$t_max,
        # Set dyn_er to 0 and overwrite during the simulation
        delta = input$t_step / input$subdiv,
        stepsize = input$t_step, # Gets overwritten anyway
        pars = list(
          a = input$a,
          omega = input$omega,
          dyn_er = sqrt(input$sigma_eta),
          meas_er = sqrt(input$sigma_epsilon),
          y_nod = input$y_nod,
          b_nod = input$b_nod,
          v_nod = input$v_nod
        ),
        start = list(
          formula(y ~ y_nod),
          formula(v ~ b_nod),
          formula(b ~ v_nod)
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
        meas_eq = list(formula(y_obs ~ y + rnorm(1, 0, meas_er)))
      )
    }
  })

  tsm_data <- eventReactive(input$sim_data, {
    if (input$dg_set_seed) {
      set.seed(input$dg_seed)
    }
    get_tsm_data(sim_tsm(gen_tsm()))
  })

  observeEvent(input$sim_data, {
    output$tsm_model_header <- renderUI(
      fluidRow(
        h3("Generated Data:")
      )
    )
  })

  observeEvent(input$sim_data, {
    output$tsm_model <- renderPlot(
      {
        plot(
          x = tsm_data()$time, y = tsm_data()$y_obs,
          axes = FALSE, xlab = "Time", ylab = "y"
        )
        lines(x = tsm_data()$time, y = tsm_data()$y)

        axis(2)

        axis(1,
          padj = 1
        )
      },
      width = 550,
      height = 400
    )
  })

  ### Model fitting ------------------------------------------------------------
  hideTab(inputId = "gen_and_switch", target = "Model Fitting")

  observeEvent(input$sim_data, {
    showTab(inputId = "gen_and_switch", target = "Model Fitting")
  })

  observeEvent(input$fit_models, {
    output$fit_ui_header <- renderUI({
      if (input$fit_lpr | input$fit_gam) {
        fluidRow(h2("Model Fitting:"))
      }
    })
  })

  tsm_data_mod_fit <- eventReactive(input$fit_models, {
    tsm_data()
  })

  ## Fit and plot local polynomial paper
  output$locpol_settings_ui <- renderUI({
    if (input$fit_lpr) {
      fluidRow(
        column(
          width = 12,
          selectInput("locpol_settings_kernel", "Kernel:",
            choices = list(
              "Gaussian" = "gau",
              "Epanechnikov" = "epa",
              "Triangular" = "tri",
              "Uniform" = "uni"
            )
          ),
          selectInput("locpol_settings_polynomial_degree", "Polynomial Degree:",
            choices = list("1", "3", "5")
          ),
          checkboxInput("lpr_summary", "Print Summary")
        )
      )
    }
  })

  locpol_inst <- eventReactive(input$fit_models, {
    if (input$fit_lpr) {
      new("method_locpol",
        method_name = "locpol",
        gen_model = input$dg_process_choice
      )
    }
  })

  locpol_fit_pap <- eventReactive(input$fit_models, {
    if (input$fit_lpr) {
      fit(locpol_inst(), tsm_data_mod_fit())
    }
  })

  observeEvent(input$fit_models, {
    output$locpol_plot_pap <- renderPlot(
      {
        if (input$fit_lpr) {
          plot(
            x = tsm_data_mod_fit()$time, y = tsm_data_mod_fit()$y_obs,
            axes = FALSE, xlab = "Time", ylab = "y"
          )
          lines(
            x = tsm_data_mod_fit()$time,
            y = tsm_data_mod_fit()$y
          )
          lines(
            x = tsm_data_mod_fit()$time,
            locpol_fit_pap()@estimate, col = "red"
          )
          lines(
            x = tsm_data_mod_fit()$time,
            locpol_fit_pap()@ci$lb, col = "red",
            lty = 2
          )
          lines(
            x = tsm_data_mod_fit()$time,
            locpol_fit_pap()@ci$ub, col = "red",
            lty = 2
          )
          axis(2)

          axis(1,
            padj = 1
          )
        }
      },
      width = 600,
      height = 400
    )
  })

  output$locpol_print_pap <- renderPrint({
    if (input$lpr_summary) print(locpol_fit_pap())
  })

  ## Fit and plot local polynomial user
  locpol_inst_usr <- eventReactive(input$fit_models, {
    if (input$fit_lpr) {
      new("method_locpol_usr",
        method_name = "locpol custom",
        gen_model = input$dg_process_choice,
        conditions = list(
          input$locpol_settings_kernel,
          input$locpol_settings_polynomial_degree
        )
      )
    }
  })

  locpol_fit_usr <- eventReactive(input$fit_models, {
    if (input$fit_lpr) {
      fit(locpol_inst_usr(), tsm_data_mod_fit())
    }
  })

  observeEvent(input$fit_models, {
    output$locpol_plot_usr <- renderPlot(
      {
        if (input$fit_lpr) {
          plot(
            x = tsm_data_mod_fit()$time, y = tsm_data_mod_fit()$y_obs,
            axes = FALSE, xlab = "Time", ylab = "y"
          )
          lines(
            x = tsm_data_mod_fit()$time,
            y = tsm_data_mod_fit()$y
          )
          lines(
            x = tsm_data_mod_fit()$time,
            locpol_fit_usr()@estimate, col = "red"
          )
          lines(
            x = tsm_data_mod_fit()$time,
            locpol_fit_usr()@ci$lb, col = "red",
            lty = 2
          )
          lines(
            x = tsm_data_mod_fit()$time,
            locpol_fit_usr()@ci$ub, col = "red",
            lty = 2
          )
          axis(2)

          axis(1,
            padj = 1
          )
        }
      },
      width = 600,
      height = 400
    )
  })

  output$locpol_print_usr <- renderPrint({
    if (input$lpr_summary) print(locpol_fit_usr())
  })

  ## Render Locpol UI
  observeEvent(input$fit_models, {
    output$locpol_ui_pop <- renderUI({
      if (input$fit_lpr) {
        fluidRow(
          h3("Local Polynomial Regression:"),
          column(
            width = 7,
            h4("LPR fit using the same configurations as in the paper:"),
            plotOutput("locpol_plot_pap", width = "auto", height = "auto"),
            verbatimTextOutput("locpol_print_pap"),
            offset = 0
          ),
          column(
            width = 5,
            h4("LPR fit using your configurations:"),
            plotOutput("locpol_plot_usr", width = "auto", height = "auto"),
            verbatimTextOutput("locpol_print_usr"),
            offset = 0
          )
        )
      }
    })
  })

  ## Fit and plot gam paper
  output$gam_settings_ui <- renderUI({
    if (input$fit_gam) {
      fluidRow(
        column(
          width = 12,
          selectInput("gam_settings_spline_basis", "Spline Basis:",
            choices = list(
              "Thin-plate" = "tp",
              "Cubic" = "cr",
              "Cyclic Cubic" = "cc",
              "P-Splines" = "ps"
            )
          ),
          numericInput("gam_settings_numbbasis",
            "Number of Basis Functions:",
            -1,
            min = -1, max = nrow(tsm_data())
          ),
          checkboxInput("gam_summary", "Print Summary")
        )
      )
    }
  })

  gam_inst <- eventReactive(input$fit_models, {
    if (input$fit_gam) {
      new("method_gam",
        method_name = "gam",
        gen_model = input$dg_process_choice
      )
    }
  })

  gam_fit_pap <- eventReactive(input$fit_models, {
    if (input$fit_gam) {
      fit(gam_inst(), tsm_data_mod_fit())
    }
  })

  observeEvent(input$fit_models, {
    output$gam_plot_pap <- renderPlot(
      {
        if (input$fit_gam) {
          plot(
            x = tsm_data_mod_fit()$time, y = tsm_data_mod_fit()$y_obs,
            axes = FALSE, xlab = "Time", ylab = "y"
          )
          lines(x = tsm_data_mod_fit()$time, y = tsm_data_mod_fit()$y)
          lines(
            x = tsm_data_mod_fit()$time,
            gam_fit_pap()@estimate, col = "red"
          )
          lines(
            x = tsm_data_mod_fit()$time,
            gam_fit_pap()@ci$lb, col = "red", lty = 2
          )
          lines(
            x = tsm_data_mod_fit()$time,
            gam_fit_pap()@ci$ub, col = "red", lty = 2
          )
          axis(2)

          axis(1,
            padj = 1
          )
        }
      },
      width = 600,
      height = 400
    )
  })

  output$gam_print_pap <- renderPrint({
    if (input$gam_summary) print(gam_fit_pap())
  })

  ## Fit and plot gam usr
  gam_inst_usr <- eventReactive(input$fit_models, {
    if (input$fit_gam) {
      new("method_gam_usr",
        method_name = "gam",
        gen_model = input$dg_process_choice,
        conditions = list(
          input$gam_settings_spline_basis,
          input$gam_settings_numbbasis
        )
      )
    }
  })

  gam_fit_usr <- eventReactive(input$fit_models, {
    if (input$fit_gam) {
      fit(gam_inst_usr(), tsm_data_mod_fit())
    }
  })

  observeEvent(input$fit_models, {
    output$gam_plot_usr <- renderPlot(
      {
        if (input$fit_gam) {
          plot(
            x = tsm_data_mod_fit()$time, y = tsm_data_mod_fit()$y_obs,
            axes = FALSE, xlab = "Time", ylab = "y"
          )
          lines(x = tsm_data_mod_fit()$time, y = tsm_data_mod_fit()$y)
          lines(
            x = tsm_data_mod_fit()$time,
            gam_fit_usr()@estimate, col = "red"
          )
          lines(
            x = tsm_data_mod_fit()$time,
            gam_fit_usr()@ci$lb, col = "red", lty = 2
          )
          lines(
            x = tsm_data_mod_fit()$time,
            gam_fit_usr()@ci$ub, col = "red", lty = 2
          )
          axis(2)

          axis(1,
            padj = 1
          )
        }
      },
      width = 600,
      height = 400
    )
  })

  output$gam_print_usr <- renderPrint({
    if (input$gam_summary) print(gam_fit_usr())
  })

  ### Render gam output UI

  observeEvent(input$fit_models, {
    output$gam_ui_pop <- renderUI({
      if (input$fit_gam) {
        fluidRow(
          h3("Generalized Additive Model:"),
          column(
            width = 7,
            h4("GAM fit using the same configurations as in the paper:"),
            plotOutput("gam_plot_pap", width = "auto", height = "auto"),
            verbatimTextOutput("gam_print_pap"),
            offset = 0
          ),
          column(
            width = 5,
            h4("GAM fit using your configurations:"),
            plotOutput("gam_plot_usr", width = "auto", height = "auto"), verbatimTextOutput("gam_print_usr"),
            offset = 0
          )
        )
      }
    })
  })
}

# Create app
shinyApp(ui = ui, server = server)

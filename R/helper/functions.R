#' ----------------------------------------------------------------------------#
#' Title: Helper Function File                                                 #
#' Author: Jan Ian Failenschmid                                                #
#' Created Date: 25-03-2024                                                    #
#' -----                                                                       #
#' Last Modified: 27-03-2024                                                   #
#' Modified By: Jan Ian Failenschmid                                           #
#' -----                                                                       #
#' Copyright (c) 2024 by Jan Ian Failenschmid                                  #
#' E-mail: J.I.Failenschmid@tilburguniveristy.edu                              #
#' -----                                                                       #
#' License: GNU General Public License v3.0 or later                           #
#' License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html           #
#' ----------------------------------------------------------------------------#



### Functions ------------------------------------------------------------------
simulate <- function(gen_model_list, method_list, time, repetitions) {
  ## Setup-simulation and initialize objects

  # Create sim object by initiating objects from the gen_model list
  sim_grid <- expand.grid(time = time, gen_model = gen_model_list)
  sim_grid$gen_model <- mapply(
    function(model, time) {
      model@time <- time
      return(model)
    },
    model = sim_grid$gen_model, time = sim_grid$time
  )

  # Replicate the rows of sim_grid repetitions time
  sim_grid <- sim_grid[rep(seq_len(nrow(sim_grid)), each = repetitions), ]

  # Add unique identifier for data files
  sim_grid$dat_id <- paste(
    sapply(
      sim_grid$gen_model,
      function(x) slot(x, "model_name")
    ),
    sim_grid$time,
    rep(seq_len(repetitions), length(time) * length(gen_model_list)),
    sep = "_"
  )

  # Initialize results grid by expanding the data ids over the elements of
  # method_list
  res_grid <- expand.grid(dat_id = sim_grid$dat_id, method = method_list)

  # Update the time and gen_model slot of the method entries accordingly
  res_grid$method <- mapply(
    function(id, method, sim_grid) {
      gen_model <- sim_grid$gen_model[[which(sim_grid$dat_id == id)]]
      method@time <- gen_model@time
      method@gen_model <- gen_model@model_name
      return(method)
    },
    id = res_grid$dat_id, method = res_grid$method,
    MoreArgs = list(sim_grid = sim_grid), SIMPLIFY = FALSE
  )

  # Sample models and fit methods until the desired amount of repetitions is
  # achieved
  while (!all(sapply(res_grid$method, function(x) x@converged) == TRUE)) {
    # Indicators for which models are converged and which data sets
    # need to be resampled
    ind_res <- !sapply(res_grid$method, function(x) x@converged)
    ind_sim <- sim_grid$dat_id %in% res_grid$dat_id[ind_res]

    ## Data simulation

    # Simulate data from each model object
    sim_grid$dat[ind_sim] <- lapply(sim_grid$gen_model[ind_sim], sim_ssm)

    ## Model fitting

    # Fit models
    res_grid$fit[ind_res] <- mapply(
      function(id, method, sim_grid) {
        fit(
          method = method,
          data = get_ssm_data(
            sim_grid$dat[[which(sim_grid$dat_id == id)]], TRUE, TRUE
          )
        )
      },
      id = res_grid$dat_id[ind_res], method = res_grid$method[ind_res],
      MoreArgs = list(sim_grid = sim_grid), SIMPLIFY = FALSE
    )

    # Copy convergence indicator into method object
    res_grid$method <- mapply(function(method, fit) {
      method@converged <- fit$converged
      return(method)
    }, method = res_grid$method, fit = res_grid$fit)
  }

  ## State inference

  # Extract state inference from fit object
  res_grid$state_inference <- mapply(function(method, fit) {
    infer_state(method, fit)
  }, method = res_grid$method, fit = res_grid$fit, SIMPLIFY = FALSE)

  ## Calculate performance measures

  # Check the proportion of time that the state falls within the confidence/
  # credible interval
  res_grid$ci_prop <- mapply(
    function(method, state_inf, id, sim_grid) {
      ci_test(
        method,
        state_inf,
        get_ssm_data(
          sim_grid$dat[[which(sim_grid$dat_id == id)]], FALSE,
          FALSE, TRUE
        )
      )
    },
    method = res_grid$method,
    state_inf = res_grid$state_inf,
    id = res_grid$dat_id,
    MoreArgs = list(sim_grid = sim_grid), SIMPLIFY = FALSE
  )

  # Calculate the mean squared error for the state inference
  res_grid$mse <- mapply(
    function(method, state_inf, id, sim_grid) {
      calc_mse(
        method,
        state_inf,
        get_ssm_data(
          sim_grid$dat[[which(sim_grid$dat_id == id)]], FALSE,
          FALSE, TRUE
        )
      )
    },
    method = res_grid$method,
    state_inf = res_grid$state_inf,
    id = res_grid$dat_id,
    MoreArgs = list(sim_grid = sim_grid), SIMPLIFY = FALSE
  )

  # Calculate GCV
  res_grid$gcv <- mapply(
    function(method, fit, state_inf, id, sim_grid) {
      calc_gcv(
        method,
        fit,
        state_inf,
        get_ssm_data(
          sim_grid$dat[[which(sim_grid$dat_id == id)]],
          TRUE, TRUE
        )
      )
    },
    method = res_grid$method, fit = res_grid$fit,
    state_inf = res_grid$state_inf, id = res_grid$dat_id,
    MoreArgs = list(sim_grid = sim_grid), SIMPLIFY = FALSE
  )

  return(list(sim_grid = sim_grid, res_grid = res_grid))
}

sim_ssm <- function(model) {
  # Function to actually simulate the data

  # Input checks and formating ----
  # If model time is not a positive integer, return error
  if (is.na(slot(model, "time")) || slot(model, "time") < 1 ||
    slot(model, "time") %% 1 != 0) {
    stop("Time needs to be a positive integer.")
  }

  # If model type is not valid, return error
  if (is.na(slot(model, "model_type")) ||
    !slot(model, "model_type") %in% c("SSM", "DE")) {
    stop("Model type needs to be either SSM or DE.")
  }

  # If delat is negative or a non-integer devisor of 1, return error
  if (slot(model, "model_type") == "DE") {
    if (is.na(slot(model, "delta")) || slot(model, "delta") <= 0) {
      stop("Delta needs to be a positive real number.")
    }
    if ((1 / slot(model, "delta")) %% 1 != 0) {
      stop("Delta needs to be a integer devisor of 1.")
    }
  }

  # If dynamic error list is empty, set dynamic error to zero for all variables
  if (length(slot(model, "dynamic_error")) == 0) {
    slot(model, "dynamic_error") <- lapply(
      lapply(slot(model, "state_eq"), "[[", 2),
      function(x) {
        as.formula(paste0(x, " ~ 0"))
      }
    )
  }

  # Match formula order between start, state_eq, and dynamic_error
  slot(model, "start") <- slot(model, "start")[match(
    sapply(slot(model, "state_eq"), "[[", 2),
    sapply(slot(model, "start"), "[[", 2)
  )]

  slot(model, "dynamic_error") <- slot(model, "dynamic_error")[match(
    sapply(slot(model, "state_eq"), "[[", 2),
    sapply(slot(model, "dynamic_error"), "[[", 2)
  )]

  # If any dynamic error is negative, return error
  if (!all(mapply(
    function(x, y) {
      start <- eval(y[[3]], slot(model, "pars"))
      names(start) <- deparse(y[[2]])
      eval(x[[3]], list2env(c(start, slot(model, "pars"))))
    },
    x = slot(model, "dynamic_error"), y = slot(model, "start")
  ) >= 0)) {
    stop("All dyamic errors need to be real non-negative values")
  }

  # Simulation ----

  # Initialize data frame to hold simulation results
  dat <- data.frame(time = seq(1, slot(model, "time")))

  # Set initial state value to specified starting value
  dat$state[[1]] <- lapply(slot(model, "start"), function(x) {
    eval(x[[3]], slot(model, "pars"))
  })
  names(dat$state[[1]]) <- lapply(slot(model, "start"), function(x) x[[2]])

  if (slot(model, "model_type") == "SSM") {
    # Apply state function iteratively to get state values starting form the
    # second time point
    for (t in seq(2, slot(model, "time"))) {
      dat$state[[t]] <- mapply(
        # Apply function to combinations of state and dynamic error expressions
        function(x, y) {
          # Evaluate dynamic error function in parameter environment
          error <- eval(
            str2lang(
              paste(deparse(y[[3]]),
                rnorm(1, 0, 1),
                sep = " * "
              )
            ),
            list2env(c(dat$state[[t - 1]], slot(model, "pars")))
          )
          # Evaluate state function in parameter environment and add error
          state_plus_one <- eval(
            str2lang(paste(deparse(x[[3]]), error,
              sep = " + "
            )),
            list2env(c(dat$state[[t - 1]], slot(model, "pars")))
          )
          return(slot(model, "boundary")(state_plus_one))
        },
        x = slot(model, "state_eq"), y = slot(model, "dynamic_error"),
        SIMPLIFY = FALSE
      )
      # Add names
      names(dat$state[[t]]) <- lapply(
        slot(model, "state_eq"),
        function(x) x[[2]]
      )
    }
  } else if (slot(model, "model_type") == "DE") {
    # Run Euler-Maruyama method starting from second time point
    # Instantiate new state list that holds the oversampled state
    state <- list(dat$state[[1]])
    for (t in seq_len(slot(model, "time") / slot(model, "delta"))[-1]) {
      state[[t]] <- mapply(
        # Apply function to combinations of state and dynamic error expressions
        function(x, y) {
          # Evaluate the Wiener (diffusion) component in the parameter env.
          wiener <- eval(
            str2lang(
              paste(deparse(y[[3]]),
                rnorm(1, 0, sqrt(slot(model, "delta"))),
                sep = " * "
              )
            ),
            list2env(c(state[[t - 1]], slot(model, "pars")))
          )
          # Evaluate the differential (drift) component in the parameter env
          # and add previous state value und Wiender component
          state_plus_one <- eval(
            str2lang(paste0(
              deparse(x[[2]]), " + ", slot(model, "delta"), "*(",
              deparse(x[[3]]), ")", " + ", wiener
            )),
            list2env(c(state[[t - 1]], slot(model, "pars")))
          )
          return(slot(model, "boundary")(state_plus_one))
        },
        x = slot(model, "state_eq"), y = slot(model, "dynamic_error"),
        SIMPLIFY = FALSE
      )
      # Add names
      names(state[[t]]) <- lapply(
        slot(model, "state_eq"),
        function(x) x[[2]]
      )
    }
    # Copy subsample at the desired time points over to data frame
    dat$state <- state[seq_len(slot(model, "time")) *
      (1 / slot(model, "delta")) - (1 / slot(model, "delta")) + 1]
  }

  # Generate observation values from state values
  dat$observation <- lapply(
    dat$state,
    function(state) {
      val <- lapply(
        slot(model, "meas_eq"),
        function(x) {
          eval(
            x[[3]],
            list2env(c(state, slot(model, "pars")))
          )
        }
      )
      names(val) <- lapply(slot(model, "meas_eq"), function(x) x[[2]])
      return(val)
    }
  )

  return(dat)
}

get_ssm_data <- function(data,
                         time = FALSE,
                         observation = FALSE,
                         state = FALSE) {
  # Function to extract data features from sim_ssm generated list

  # Validation
  if (!time && !observation && !state) {
    cat("Select at least one data feature to be included. \n")
  }

  # Select data featrues
  res <- list()

  if (time) res$time <- data$time
  if (observation) {
    res$observation <- do.call(
      rbind.data.frame,
      lapply(data$observation, rbind)
    )
  }
  if (state) res$state <- do.call(rbind.data.frame, lapply(data$state, rbind))

  nam <- unlist(sapply(res, names)) # Extract names
  if (time) nam <- c("time", nam) # Add time name
  res <- do.call(cbind, res) # Reformat data
  names(res) <- nam # Fix names

  res <- as.data.frame(sapply(res, unlist)) # Get rid of lists within columns

  return(res)
}

ci_test <- function(method, state_inf, data) {
  ci_incl <- mapply(function(y, ub, lb) y < ub & y > lb,
    y = data$y, ub = state_inf$ub, lb = state_inf$lb
  )

  smooth <- mean(ci_incl)
  return(list(ci_prop = smooth, ci_incl = ci_incl))
}

calc_mse <- function(method, state_inf, data) {
  squared_err <- (state_inf$y_hat - data$y)^2
  smooth <- mean(squared_err)
  return(smooth)
}

plot_state_inference <- function(sim, row = 1, observation, state, estimate) {
  state_inf <- sim$res_grid$state_inf[[row]]
  data <- sim$sim_grid$dat[sim$sim_grid$dat_id == sim$res_grid$dat_id[[row]]]
  data <- get_ssm_data(data[[1]], TRUE, TRUE, TRUE)


  plot(x = data$time, data[[observation]])
  lines(x = data$time, y = data[[state]])
  lines(x = data$time, y = state_inf[[estimate]], col = "red")
  lines(x = data$time, y = state_inf$lb, col = "red", lty = 2)
  lines(x = data$time, y = state_inf$ub, col = "red", lty = 2)
}

extract_results <- function(sim, col_names) {
  res_grid <- sim$res_grid
  results <- lapply(res_grid$method, function(x) {
    list(
      method = x@method_name,
      model = x@gen_model,
      time = x@time
    )
  })

  results <- do.call(rbind.data.frame, results)

  results$method <- as.factor(results$method)
  results$model <- as.factor(results$model)
  results$time <- as.factor(results$time)

  for (col in col_names) {
    if (class(res_grid[[col]][[1]]) == "list") {
      res <- list()
      for (i in seq_len(length(res_grid[[col]][[1]]))) {
        res[[i]] <- lapply(res_grid[[col]], function(x, i) x[[i]], i = i)
      }
      res <- do.call(cbind, res)
      colnames(res) <- paste(col, names(res_grid[[col]][[1]]), sep = "_")
      results <- cbind(results, res)
    } else {
      results <- cbind(results, unlist(res_grid[[col]]))
      names(results)[ncol(results)] <- col
    }
  }

  return(results)
}

# plot_mlts <- function(
#     sim_grid,
#     gen_model,
#     method,
#     time,
#     quant = c(.80, .95, .975),
#     indv = FALSE,
#     poly = TRUE,
#     var = "smooth",
#     state_name = NULL) {
#   # Function for plotting multilevel scalar time series as an output of the simulation
#   # function by generative model, analysis method, variable, and state name.

#   # For var prepare data for plotting
#   if (var == "smooth_error") {
#     # Calculate smoothing error by extracting the smoothed value and the observation
#     # and differencing
#     sim_grid$smooth_error <- mapply(
#       function(dat, smooth) {
#         dat <- getSSMdata(dat, what = "observation", time = TRUE)
#         return(data.frame(
#           time = dat$time,
#           error = smooth[, names(smooth) != "time"] -
#             dat[, names(smooth) != "time"]
#         ))
#       },
#       dat = sim_grid$dat, sim_grid$smooth, SIMPLIFY = FALSE
#     )
#   } else if (var == "state") {
#     # Extract state by state name from the dat object in sim_grid
#     sim_grid$state <- lapply(
#       sim_grid$dat,
#       function(dat, what, time, state_name) {
#         getSSMdata(dat, what = what, time = time)[, c("time", state_name)]
#       },
#       what = "state", time = TRUE, state_name = state_name
#     )
#   } else if (var == "observation") {
#     # Extract observation value from the dat object in sim grid
#     sim_grid$observation <- lapply(sim_grid$dat, getSSMdata,
#       what = "observation",
#       time = TRUE
#     )
#   }

#   # subset sim_grid by gen_model, method, time, and var
#   dat_list <- sim_grid[sim_grid$gen_model == gen_model &
#     sim_grid$method == method &
#     sim_grid$time == time, var]
#   dat <- do.call(rbind.data.frame, dat_list)

#   # Initialize plot_dat data frame with time column
#   plot_dat <- data.frame(time = unique(dat$time))

#   # Add mean column by averaging var at each time point
#   plot_dat$mean <- sapply(plot_dat$time,
#     function(time, dat) {
#       mean(dat[dat$time == time, names(dat) != "time"])
#     },
#     dat = dat
#   )

#   # Add columns for the bounds of the polynomials depending on the 3 qunatile sizes
#   plot_dat$lb1 <- sapply(plot_dat$time,
#     function(time, dat, quant) {
#       quantile(dat[dat$time == time, names(dat) != "time"], quant)
#     },
#     dat = dat, quant = 1 - quant[1]
#   )

#   plot_dat$ub1 <- sapply(plot_dat$time,
#     function(time, dat, quant) {
#       quantile(dat[dat$time == time, names(dat) != "time"], quant)
#     },
#     dat = dat, quant = quant[1]
#   )

#   plot_dat$lb2 <- sapply(plot_dat$time,
#     function(time, dat, quant) {
#       quantile(dat[dat$time == time, names(dat) != "time"], quant)
#     },
#     dat = dat, quant = 1 - quant[2]
#   )

#   plot_dat$ub2 <- sapply(plot_dat$time,
#     function(time, dat, quant) {
#       quantile(dat[dat$time == time, names(dat) != "time"], quant)
#     },
#     dat = dat, quant = quant[2]
#   )

#   plot_dat$lb3 <- sapply(plot_dat$time,
#     function(time, dat, quant) {
#       quantile(dat[dat$time == time, names(dat) != "time"], quant)
#     },
#     dat = dat, quant = 1 - quant[3]
#   )

#   plot_dat$ub3 <- sapply(plot_dat$time,
#     function(time, dat, quant) {
#       quantile(dat[dat$time == time, names(dat) != "time"], quant)
#     },
#     dat = dat, quant = quant[3]
#   )

#   if ("forcast" %in% names(sim_grid)) {
#     # Repeat everything for forcast if forcast is included

#     # For var prepare data for plotting
#     if (var == "smooth_error") {
#       # Calculate smoothing error by extracting the smoothed value and the observation
#       # and differencing
#       sim_grid$forc_smooth_error <- mapply(
#         function(forc_dat, forcast) {
#           forc_dat <- getSSMdata(forc_dat, what = "observation", time = TRUE)
#           return(data.frame(
#             time = forc_dat$time,
#             error = forcast[, names(forcast) != "time"] -
#               forc_dat[, names(forcast) != "time"]
#           ))
#         },
#         forc_dat = sim_grid$forc_dat, sim_grid$forcast, SIMPLIFY = FALSE
#       )
#     } else if (var == "state") {
#       # Extract state by state name from the dat object in sim_grid
#       sim_grid$forc_state <- lapply(
#         sim_grid$forc_dat,
#         function(forc_dat, what, time, state_name) {
#           getSSMdata(
#             forc_dat,
#             what = what,
#             time = time
#           )[, c("time", state_name)]
#         },
#         what = "state", time = TRUE, state_name = state_name
#       )
#     } else if (var == "observation") {
#       # Extract observation value from the dat object in sim grid
#       sim_grid$forc_observation <- lapply(sim_grid$forc_dat, getSSMdata,
#         what = "observation",
#         time = TRUE
#       )
#     }

#     if (var == "smooth") {
#       # subset sim_grid by gen_model, method, time, and forc_var
#       forc_dat_list <- sim_grid[sim_grid$gen_model == gen_model &
#         sim_grid$method == method &
#         sim_grid$time == time, "forcast"]
#     } else {
#       # subset sim_grid by gen_model, method, time, and forc_var
#       forc_dat_list <- sim_grid[sim_grid$gen_model == gen_model &
#         sim_grid$method == method &
#         sim_grid$time == time, paste0("forc_", var)]
#     }

#     forc_dat <- do.call(rbind.data.frame, forc_dat_list)

#     # Initialize plot_forc_dat data frame with time column
#     plot_forc_dat <- data.frame(time = unique(forc_dat$time))

#     # Add mean column by averaging var at each time point
#     plot_forc_dat$mean <- sapply(plot_forc_dat$time,
#       function(time, forc_dat) {
#         mean(forc_dat[forc_dat$time == time, names(forc_dat) != "time"])
#       },
#       forc_dat = forc_dat
#     )

#     # Add columns for the bounds of the polynomials depending on the 3 qunatile sizes
#     plot_forc_dat$lb1 <- sapply(
#       plot_forc_dat$time,
#       function(time, forc_dat, quant) {
#         quantile(
#           forc_dat[
#             forc_dat$time == time,
#             names(forc_dat) != "time"
#           ],
#           quant
#         )
#       },
#       forc_dat = forc_dat, quant = 1 - quant[1]
#     )

#     plot_forc_dat$ub1 <- sapply(
#       plot_forc_dat$time,
#       function(time, forc_dat, quant) {
#         quantile(
#           forc_dat[
#             forc_dat$time == time,
#             names(forc_dat) != "time"
#           ],
#           quant
#         )
#       },
#       forc_dat = forc_dat, quant = quant[1]
#     )

#     plot_forc_dat$lb2 <- sapply(
#       plot_forc_dat$time,
#       function(time, forc_dat, quant) {
#         quantile(
#           forc_dat[
#             forc_dat$time == time,
#             names(forc_dat) != "time"
#           ],
#           quant
#         )
#       },
#       forc_dat = forc_dat, quant = 1 - quant[2]
#     )

#     plot_forc_dat$ub2 <- sapply(
#       plot_forc_dat$time,
#       function(time, forc_dat, quant) {
#         quantile(
#           forc_dat[
#             forc_dat$time == time,
#             names(forc_dat) != "time"
#           ],
#           quant
#         )
#       },
#       forc_dat = forc_dat, quant = quant[2]
#     )

#     plot_forc_dat$lb3 <- sapply(
#       plot_forc_dat$time,
#       function(time, forc_dat, quant) {
#         quantile(
#           forc_dat[
#             forc_dat$time == time,
#             names(forc_dat) != "time"
#           ],
#           quant
#         )
#       },
#       forc_dat = forc_dat, quant = 1 - quant[3]
#     )

#     plot_forc_dat$ub3 <- sapply(
#       plot_forc_dat$time,
#       function(time, forc_dat, quant) {
#         quantile(
#           forc_dat[
#             forc_dat$time == time,
#             names(forc_dat) != "time"
#           ],
#           quant
#         )
#       },
#       forc_dat = forc_dat, quant = quant[3]
#     )
#   }

#   # Generate plot main title based on selected var
#   if (var == "smooth") {
#     main_str <- paste0(
#       "Smoothed observations: ", gen_model,
#       " ", method, " ", time
#     )
#   } else if (var == "smooth_error") {
#     main_str <- paste0("Smoothing error: ", gen_model, " ", method, " ", time)
#   } else if (var == "observation") {
#     main_str <- paste0("Observations: ", gen_model, " ", method, " ", time)
#   } else if (var == "state") {
#     main_str <- paste0("State: ", gen_model, " ", method, " ", time)
#   }

#   if ("forcast" %in% names(sim_grid)) {
#     # Initialize plot with the correct labels and dimensions
#     plot(1,
#       type = "n", xlab = "Time", ylab = "Observation",
#       main = main_str,
#       xlim = c(
#         min(c(plot_dat$time, plot_forc_dat$time)),
#         max(c(plot_dat$time, plot_forc_dat$time))
#       ),
#       ylim = c(
#         min(c(plot_dat$lb3, plot_forc_dat$lb3)),
#         max(c(plot_dat$ub3, plot_forc_dat$ub3))
#       )
#     )
#   } else {
#     # Initialize plot with the correct labels and dimensions
#     plot(1,
#       type = "n", xlab = "Time", ylab = "Observation",
#       main = main_str,
#       xlim = c(min(plot_dat$time), max(plot_dat$time)),
#       ylim = c(min(plot_dat$lb3), max(plot_dat$ub3))
#     )
#   }

#   # If polygon should be plotted, add polygons
#   if (poly) {
#     # Add outer ploygon
#     polygon(
#       x = c(plot_dat$time, rev(plot_dat$time)),
#       y = c(plot_dat$lb3, rev(plot_dat$ub3)),
#       col = "lightblue"
#     )

#     # Add middle polygon
#     polygon(
#       x = c(plot_dat$time, rev(plot_dat$time)),
#       y = c(plot_dat$lb2, rev(plot_dat$ub2)),
#       col = "blue"
#     )

#     # Add inner polygon
#     polygon(
#       x = c(plot_dat$time, rev(plot_dat$time)),
#       y = c(plot_dat$lb1, rev(plot_dat$ub1)),
#       col = "darkblue"
#     )

#     # If polygon should be plotted, add polygons
#     if ("forcast" %in% names(sim_grid)) {
#       # Add outer ploygon
#       polygon(
#         x = c(plot_forc_dat$time, rev(plot_forc_dat$time)),
#         y = c(plot_forc_dat$lb3, rev(plot_forc_dat$ub3)),
#         col = "lightgreen"
#       )

#       # Add middle polygon
#       polygon(
#         x = c(plot_forc_dat$time, rev(plot_forc_dat$time)),
#         y = c(plot_forc_dat$lb2, rev(plot_forc_dat$ub2)),
#         col = "green"
#       )

#       # Add inner polygon
#       polygon(
#         x = c(plot_forc_dat$time, rev(plot_forc_dat$time)),
#         y = c(plot_forc_dat$lb1, rev(plot_forc_dat$ub1)),
#         col = "darkgreen"
#       )
#     }
#   }

#   # Add line for the mean data
#   lines(x = plot_dat$time, y = plot_dat$mean)

#   # Add points for the mean data
#   points(x = plot_dat$time, y = plot_dat$mean)

#   if ("forcast" %in% names(sim_grid)) {
#     # Add line for the mean data
#     lines(x = plot_forc_dat$time, y = plot_forc_dat$mean)

#     # Add points for the mean data
#     points(x = plot_forc_dat$time, y = plot_forc_dat$mean)
#   }

#   # Plot individual observations as lines, if desired
#   if (indv) {
#     lapply(dat_list, function(x) {
#       lines(x$time, x[, names(x) != "time"], col = rgb(0.5, 0.5, 0.5, 0.3))
#     })

#     if ("forcast" %in% names(sim_grid)) {
#       lapply(forc_dat_list, function(x) {
#         lines(x$time, x[, names(x) != "time"], col = rgb(0.5, 0.5, 0.5, 0.3))
#       })
#     }
#   }
# }

gp_post_pred <- function(
    fit, time, obs, state, f_name, draws = 200,
    alpha = 1, lwd = 1, rect = FALSE,
    quant = c(.60, .80, .95)) {
  # Function to do posterior predictive plotting of gaussian processes

  f_draws <- as.data.frame(fit$draws(f_name, format = "draws_df"))
  f_draws <- f_draws[, seq(1, ncol(f_draws) - 3)]

  if (!rect) {
    plot(x = time, y = obs)

    lapply(sample(seq(1, nrow(f_draws)), draws),
      function(x, f_draws, time, lwd, alpha) {
        lines(
          x = time, y = f_draws[x, time],
          col = rgb(red = 1, green = 0, blue = 0, alpha = alpha),
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

gp_posterior_3d <- function(
    fit, f_name, draws = 200, grid_size = 225,
    state = NULL, obs = NULL, smooth = NULL) {
  require(plotly)
  require(MASS)
  require(tidyr)

  f_draws <- as.data.frame(fit$draws("f", format = "draws_df"))
  f_draws <- f_draws[
    sample(seq_len(nrow(f_draws)), draws),
    seq(1, ncol(f_draws) - 3)
  ]
  time <- seq_len(ncol(f_draws))
  colnames(f_draws) <- time
  f_draws <- gather(f_draws, key = "time", value = "state")
  f_draws$time <- as.numeric(f_draws$time)

  kd <- kde2d(f_draws$time, f_draws[["state"]], n = 225)

  fig <- plot_ly(x = kd$x, y = kd$y, z = t(kd$z)) %>%
    add_surface()

  if (!is.null(state)) {
    z <- rep(max(kd$z) * 1.1, length(state))
    fig <- add_trace(fig,
      y = ~state,
      x = ~time,
      z = ~z,
      type = "scatter3d", mode = "lines",
      line = list(size = 5, color = "red")
    )
  }

  if (!is.null(smooth)) {
    z <- rep(max(kd$z) * 1.1, length(smooth))
    fig <- add_trace(fig,
      y = ~smooth,
      x = ~time,
      z = ~z,
      type = "scatter3d", mode = "lines",
      line = list(size = 5, color = "#00eeff")
    )
  }

  if (!is.null(obs)) {
    z <- rep(max(kd$z) * 1.1, length(state))
    fig <- add_markers(fig,
      y = ~obs,
      x = ~time,
      z = ~z,
      type = "scatter3d",
      marker = list(size = 3, color = "orange")
    )
  }

  fig <- layout(fig,
    showlegend = FALSE,
    title = "Posterior Gaussian Process Distribution",
    scene = list(
      xaxis = list(title = "Time"),
      yaxis = list(title = "State"),
      zaxis = list(title = "Probability Density"),
      camera = list(eye = list(x = 1.25, y = 1.25, z = 1.25))
    )
  )

  fig
}

make_exemplar_plot <- function(gen_model, main, xlab, ylab, cex) {
  dat_list <- sim_ssm(gen_model)
  dat <- get_ssm_data(dat_list, TRUE, TRUE, TRUE)

  plot(dat$time, dat$y_obs,
    main = main, xlab = xlab, ylab = ylab,
    cex.lab = cex, cex.axis = cex, cex.main = cex, cex.sub = cex
  )
  lines(dat$time, dat$y)
}

#' ----------------------------------------------------------------------------#
#' Title: Classes and Generic Methods                                          #
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



### Define classes, generics and methods ---------------------------------------

## gen_model

# Classes

setClass(
  "gen_model",
  slots = list(
    model_name = "character",
    model_type = "character",
    delta = "numeric",
    dynamic_error = "list",
    time = "numeric",
    pars = "list",
    start = "list",
    state_eq = "list",
    meas_eq = "list",
    boundary = "function"
  ),
  prototype = list(
    time = NA_integer_,
    model_type = NA_character_,
    delta = NA_real_,
    dynamic_error = list(),
    boundary = function(x) ifelse(is.finite(x), x, NA)
  )
)

# Classes

#' Method class is a convenience class to hold the structure for all method
#' subclasses.
#' It should never be initialized, instead initialize any subclass to have
#' specified methods.
setClass(
  "method",
  slots = list(
    method_name = "character",
    gen_model = "character",
    time = "numeric",
    pars = "list",
    converged = "logical"
  ),
  prototype = list(
    method_name = NA_character_,
    gen_model = NA_character_,
    time = NA_integer_,
    pars = list(),
    converged = FALSE
  )
)

setClass(
  "gam",
  contains = "method"
)

setClass(
  "locpol",
  contains = "method"
)

setClass(
  "gp",
  contains = "method"
)

setClass(
  "ssm",
  contains = "method"
)

setClass(
  "locpol2",
  contains = "method"
)

setClass(
  "locpol3",
  contains = "method"
)

# Generics
setGeneric("fit", function(method, data) standardGeneric("fit"),
  signature = "method"
)

setGeneric("infer_state", function(method, fit) standardGeneric("infer_state"),
  signature = "method"
)

setGeneric("calc_gcv",
  function(method, fit, state_inf, data) standardGeneric("calc_gcv"),
  signature = "method"
)

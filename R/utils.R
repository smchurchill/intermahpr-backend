## intermahpr - R package backend for the intermahp shiny app
## Copyright (C) 2018 Canadian Institute for Substance Use Research

#### Data carpentry ------------------------------------------------------------

#' Clean data
#'
#' Performs universal cleaning operations of lowering variable name case,
#' ensuring necessary variables, imputing missing values when possible, and
#' reconverting types after imputation.
#'
#' @param .data a data_frame object
#' @param expected character vector of expected variables
#'
#' @return The cleaned dataset
#'
#' @export
clean <- function(.data, expected) {
  .data %<>%
    lowerVars() %>%
    checkVars(expected) %>%
    imputeMissing()

  suppressMessages(readr::type_convert(.data))
}

#' Converts data variable names to lower case and returns data
#'
#' Lowers the case of all variable names.  For compatibility with SAS version.
#'
#' @param .data a data_frame object
#'
#' @return The dataset with variable names in lower case
#'
lowerVars <- function(.data) {
  names(.data) <- stringr::str_to_lower(names(.data))
  .data
}

#' Checks variable names
#'
#' Checks against the given list, errors if a variable is missing and returns
#' only the specified variables.
#'
#' @inheritParams clean
#'
checkVars <- function(.data, expected) {
  missing <- expected[!(expected %in% names(.data))]
  if(length(missing) > 0) {

    ## This logic is getting a bit convoluted.
    ##
    ## The idea behind this is to make it easier to adapt intermahp to future
    ## research in distributions that describe group drinking patterns.
    ## Currently, if there are no gamma_constant or gamma_stderror variables,
    ## but the supplied genders are "Female" and "Male", intermahpr will supply
    ## relationships for the assumed gamma distribution.
    if(length(missing) == 2 && setequal(missing, c("gamma_constant", "gamma_stderror")) && setequal(.data$gender, c("Male", "Female"))){
      gc = list("Female" = 1.258, "Male" = 1.171) # Means, Kehoe et al. (2012)
      gs = list("Female" = (1.293-1.223) / 2 / 1.96,
                "Male" = (1.197 - 1.144) / 2 / 1.96) # Bounds on 95% CI, Kehoe et al. (2012)
      .data %<>% mutate(
        gamma_constant = map_dbl(gender, ~`[[`(gc, .x)),
        gamma_stderror = map_dbl(gender, ~`[[`(gs, .x))
      )
    } else if (length(missing) == 1 && setequal(missing, c("gamma_constant")) && setequal(.data$gender, c("Male", "Female"))) {
      gc = list("Female" = 1.258, "Male" = 1.171) # Means, Kehoe et al. (2012)
      .data %<>% mutate(
        gamma_constant = map_dbl(gender, ~`[[`(gc, .x))
      )
    } else if (length(missing) == 1 && setequal(missing, c("gamma_stderror")) && setequal(.data$gender, c("Male", "Female"))) {
      gs = list("Female" = (1.293-1.223) / 2 / 1.96,
                "Male" = (1.197 - 1.144) / 2 / 1.96) # Bounds on 95% CI, Kehoe et al. (2012)
      .data %<>% mutate(
        gamma_stderror = map_dbl(gender, ~`[[`(gs, .x))
      )
    } else {
      missing <- setdiff(missing,  c("gamma_constant", "gamma_stderror"))
      msg <- c("A supplied table was missing necessary variables:\n", paste0("\t", missing, "\n"))
      message(msg)
      stop(msg)
    }
  }
  .data[expected]
}

#' Provides a value to impute with
#'
#' @param var a variable name (as a character vector).
#'
#' @return An error if the variable is not suitable or imputation or the
#' imputation value if the variable is
#'
#'
imputeWith <- function(var) {
  if(grepl("[0-9]", var)) {var = "beta"}
  switch(
    var,
    rr_fd = 1,
    bingef = 1,
    beta = 0,
    count = 0
  )
}

#' Impute missing values over entire dataset
imputeMissing <- function(.data) {
  for(var in names(.data)) {
    .data[var][isMissing(.data[var])] <- imputeWith(var)
  }
  .data
}

#' Define missing data
#'
isMissing <- function(obs) {
  is.na(obs) | is.null(obs) | obs == "."
}

#' Get variables expected to be in the given object type
#'
#' @export
getExpectedVars <- function(...) {
  unlist(map(list(...), getExpectedVars_))
}

#' Get variables expected to be in the given object type
#'
getExpectedVars_ <- function(.obj_type) {
  switch(
    .obj_type,
    rr = c(
      "im",
      "condition",
      "gender",
      "outcome",
      "rr_fd",
      "bingef",
      "form",
      "attributability",
      paste0("b", 1:16)
    ),
    pc = c(
      "region",
      "year",
      "gender",
      "gamma_constant",
      "gamma_stderror",
      "age_group",
      "population",
      "pcc_litres_year",
      "correction_factor",
      "relative_consumption",
      "p_la",
      "p_fd",
      "p_cd",
      "p_bd"
    ),
    pc_display = c(
      "region",
      "year",
      "gender",
      "age_group",
      "population",
      "pcc_among_popn",
      "gamma_shape",
      "gamma_scale",
      "nc",
      "p_la",
      "p_fd",
      "p_cd",
      "p_bd"
    ),
    mm = c(
      "im",
      "region",
      "year",
      "gender",
      "age_group",
      "outcome",
      "count"
    ),
    model = c(
      "region",
      "year",
      "gender",
      "age_group",
      "im",
      "condition",
      "outcome",
      "incidence",
      "attributability",
      "current_fraction_factory",
      "former_fraction_factory"
    ),
    scenario = c(
      "region",
      "year",
      "gender",
      "age_group",
      "im",
      "condition",
      "outcome",
      "attributability",
      "current_fraction",
      "former_fraction"
    ),
    constants = c(
      "bb",
      "lb",
      "ub"
    ),
    distill_by = c(
      "region",
      "year",
      "gender",
      "age_group",
      "im",
      "condition",
      "outcome",
      "attributability"
    )
  )
}

#### Factories -----------------------------------------------------------------

#' Factory for integrators
makeIntegrator <- function(f, lb, ub) {
  integrate_up_to <- function(to) {
    if(to <= lb) to = lb
    if(to >= ub) to = ub
    0.1 * integrate(f = function(x) 10 * f(x), lower = lb, upper = to, abs.tol = 0L, subdivisions = 250L)$value
  }

  function(x) {
    vapply(x, integrate_up_to, 0)
  }
}

#' Factory for Pointwise Function Products
#'@description factory that produces the product of a pair of functions, where
#'the product used is pointwise multiplication
#'
#'@param f,g function that takes a single argument and produces a value that is
#'a valid argument for the `*` function
#'
makeProduct <- function(f, g) {
  function(x) f(x) * g(x)
}

#' Binary operator for Pointwise Function Products
#'
#'@description binary operator for product_factory
#'@param f,g function that takes a single argument and produces a value that is
#'a valid argument for the `*` function
#'

`%prod%` <- function(f,g) makeProduct(f,g)

#### Imports -------------------------------------------------------------------

#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr
#' mutate filter group_by ungroup inner_join select bind_rows left_join
#' @importFrom purrr map pmap map2_dbl map_dbl map2
#' @keywords internal
foo <- function() {return("bar")}

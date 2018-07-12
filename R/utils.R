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
  .data %>%
    lowerVars() %>%
    checkVars(expected) %>%
    imputeMissing() %>%
    readr::type_convert()
}

#' Converts data variable names to lower case and returns data
#'
#' Lowers the case of all variable names.  For compatibility with SAS version.
#'
#' @param .data a data_frame object
#'
#' @return The dataset with variable names in lower case
#'
#' @export
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
#' @export
checkVars <- function(.data, expected) {
  missing <- expected[!(expected %in% names(.data))]
  if(length(missing) > 0) {
    message <- paste(
      "The following variables were expected: ",
      paste(expected, collapse = ", "),
      "\n",
      "The following variables were supplied: ",
      paste(names(.data), collapse = ", ")
    )
    stop(message)
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
#' @export
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
#' @export
imputeMissing <- function(.data) {
  for(var in names(.data)) {
    .data[var][isMissing(.data[var])] <- imputeWith(var)
  }
  .data
}

#' Define missing data
#'
#' @export
isMissing <- function(obs) {
  is.na(obs) | is.null(obs) | obs == "."
}

#' Get variables expected to be in the given object type
#'
#'
#' @export
getExpectedVars <- function(...) {
  unlist(map(list(...), getExpectedVars_))
}

#' Get variables expected to be in the given object type
#'
#'
#' @export
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
      "pcc_among_drinkers",
      "gamma_shape",
      "gamma_scale",
      "nc",
      "p_la",
      "p_fd",
      "p_cd",
      "p_bd"
    ),
    dh = c(
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
#' @export
makeIntegrator <- function(f, lb, ub) {
  integrate_up_to <- function(to) {
    if(to <= lb) to = lb
    if(to >= ub) to = ub
    integrate(f = f, lower = lb, upper = to)$value
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
#' @export
makeProduct <- function(f, g) {
  function(x) f(x) * g(x)
}

#' Binary operator for Pointwise Function Products
#'
#'@description binary operator for product_factory
#'
#'@describeIn makeProduct
#'
#'@inheritParams makeProduct
#' @export
`%prod%` <- function(f,g) makeProduct(f,g)

#### Imports -------------------------------------------------------------------

#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr
#' mutate filter group_by ungroup inner_join select bind_rows left_join
#' @importFrom purrr map pmap map2_dbl map_dbl map2
#'
foo <- function() {return("bar")}

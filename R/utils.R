#### Functions -----------------------------------------------------------------

clean <- function(.data, expected, ...) {
  .data %>%
    vars_to_lower() %>%
    check_vars(expected) %>%
    impute_missing() %>%
    readr::type_convert()
}

makeIntegrater <- function(f, lb) {
  integrate_up_to(to) {
    if(to < lb) return(0)
    integrate(f = f, lower = lb, upper = to)$value
  }

  function(x) {
    vapply(x, integrate_up_to, 0)
  }
}

#' Converts data variable names to lower case and returns data
#'

lowerVars <- function(.data) {
  names(.data) <- stringr::str_to_lower(names(.data))
  .data
}

#' Checks variable names against the given list, errors if a variable is missing
#' and returns only the specified variables.
#'
#'

checkVars <- function(.data, expected) {
  missing <- expected[!(expected %in% names(.data))]
  if(length(missing) > 0) {
    stop(
      "The following variables were expected: ",
      paste(expected, collapse = ", "),
      "The following variables were missing: ",
      paste(missing, collapse = ", ")
    )
  }
  .data[expected]
}


#'Pointwise Function Product Factory
#'@description factory that produces the product of a pair of functions, where
#'the product used is pointwise multiplication
#'
#'@param f,g function that takes a single argument and produces a value that is
#'a valid argument for the `*` function
#'

makeProduct <- function(f, g) {function(x) f(x) * g(x)}

#'Pointwise Function Product Factory as Binary Operator
#'@description binary operator for product_factory
#'
#'@inheritParams product_factory

`%prod%` <- function(f,g) makeProduct(f,g)

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' What value is imputed when a variable is missing?
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
    RR_FD = 1,
    BINGEF = 1,
    BETA = 0,
    COUNT = 0
  )
}

imputeMissing <- function(.data) {
  for(var in names(.data)) {
    .data[var][is.missing(.data[var])] <- impute_with(var)
  }
}
#' Definition of missing data for the purposes of intermahpr
#'

isMissing <- function(obs) {
  is.na(obs) | is.null(obs) | obs == "."
}

#' Dummy function for allocating memory
#'

zero <- list(fn = function(...) 0)


#### Functions -----------------------------------------------------------------

#' Converts data variable names to lower case and returns data
#'
vars_to_lower <- function(.data) {
  names(.data) <- stringr::str_to_lower(names(.data))
  .data
}

#' Checks variable names against the given list, errors if a variable is missing
#' and returns only the specified variables.
#'
#'

check_vars <- function(.data, expected) {
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

product_factory <- function(f, g) {
  function(x) f(x) * g(x)
}

#'Pointwise Function Product Factory as Binary Operator
#'@description binary operator for product_factory
#'
#'@inheritParams product_factory

`%prod%` <- function(f,g) {
  product_factory(f,g)
}

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' What value is imputed when a variable is missing?
#'
#' @param var is a string -- a variable name.
#'
#' @return A string or integer that is appropriate for imputation.
#'
#'

impute_with <- function(var) {
  if(grepl("[0-9]", var)) {var = "beta"}
  switch(
    var,
    RR_FD = 1,
    BINGEF = 1,
    BETA = 0,
    COUNT = 0
  )
}

#' Definition of missing data for the purposes of intermahpr
#'

is.missing <- function(obs) {
  is.na(obs) | is.null(obs) | obs == "."
}

#' Dummy function for allocating memory
#'
#'@param ... Accepts any input
#'
#'@return Always returns 0

zero <- list(fn = function(...) 0)


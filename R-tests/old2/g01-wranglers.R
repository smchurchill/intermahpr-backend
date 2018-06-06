##### g01-wranglers ############################################################
##
## Data wrangling functions called from functions in s01-wrangle
##
##
##
##

#### Handle missing data -------------------------------------------------------

#' What value is imputed when a variable is missing?
#'
#' In order to normalize code in future methods, missing data is removed either
#' through imputation when possible, or by throwing an error when the data is
#' necessary.  impute_with is called when a missing value \emph{can} be imputed.
#'
#'
#' @param var is a string -- a variable name.
#'
#' @return A string or integer that is appropriate for imputation.
#'
#'

impute_with <- function(var) {
  if(grepl("B[0-9]", var)) {var = "BETA"}
  switch(
    var,
    IM = "(0).(0)",
    CONDITION = "Unspecified",
    REGION = "Unspecified",
    YEAR = "Unspecified",
    AGE_GROUP = "Unspecified",
    OUTCOME = "Combined",
    FUNCTION = "FP",
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

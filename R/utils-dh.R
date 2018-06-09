#### Deaths/Hosps Specific Data Carpentry --------------------------------------

#' Prepare Death/Hospitalization Data
#'

prepareDH <- function(.data) {
  cleanDH(.data)
}

#' Clean Death/Hospitalization Data
#'

cleanDH <- function(.data) {
  clean(.data, dh_vars)
}

#' List of variables expected to be in a DH sheet
#'

dh_vars <- c(
  "im",
  "region",
  "year",
  "gender",
  "age_group",
  "outcome",
  "count"
)



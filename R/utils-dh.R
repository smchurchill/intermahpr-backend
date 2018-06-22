#### Deaths/Hosps Specific Data Carpentry --------------------------------------

#' Prepare Death/Hospitalization Data
#'
#' @export
prepareDH <- function(.data) {
  clean(.data, getExpectedVars("dh"))
}

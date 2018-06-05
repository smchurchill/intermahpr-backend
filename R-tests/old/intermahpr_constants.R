#' Get and Set InterMAHPR environment/state variables
#'

intermahpr_constants <- new.env(parent = emptyenv())

#' Get and Set upper and lower bounds of integration
#'


get_upper_bound <- function() {
  intermahpr_constants$upper_bound
}

get_lower_bound <- function() {
  intermahpr_constants$lower_bound
}

set_upper_bound <- function(value) {
  old <- intermahpr_constants$upper_bound
  intermahpr_constants$upper_bound <- value
  invisible(old)
}

set_lower_bound <- function(value) {
  old <- intermahpr_constants$lower_bound
  intermahpr_constants$lower_bound <- value
  invisible(old)
}

set_upper_bound(250)
set_lower_bound(0.03)

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

#' Get and Set drinking group bounds
#'


get_female_light_moderate_barrier <- function() {
  intermahpr_constants$flm
}

get_female_moderate_heavy_barrier <- function() {
  intermahpr_constants$fmh
}

get_female_binge_barrier <- function() {
  intermahpr_constants$fbb
}

set_female_light_moderate_barrier <- function(value) {
  old <- intermahpr_constants$flm
  intermahpr_constants$flm <- value
  invisible(old)
}

set_female_moderate_heavy_barrier <- function(value) {
  old <- intermahpr_constants$fmh
  intermahpr_constants$fmh <- value
  invisible(old)
}

set_female_binge_barrier <- function(value) {
  old <- intermahpr_constants$fbb
  intermahpr_constants$fbb <- value
  invisible(old)
}

set_female_light_moderate_barrier(10)
set_female_moderate_heavy_barrier(20)
set_female_binge_barrier(53.8)

get_male_light_moderate_barrier <- function() {
  intermahpr_constants$mlm
}

get_male_moderate_heavy_barrier <- function() {
  intermahpr_constants$mmh
}

get_male_binge_barrier <- function() {
  intermahpr_constants$mbb
}

set_male_light_moderate_barrier <- function(value) {
  old <- intermahpr_constants$mlm
  intermahpr_constants$mlm <- value
  invisible(old)
}

set_male_moderate_heavy_barrier <- function(value) {
  old <- intermahpr_constants$mmh
  intermahpr_constants$mmh <- value
  invisible(old)
}

set_male_binge_barrier <- function(value) {
  old <- intermahpr_constants$mbb
  intermahpr_constants$mbb <- value
  invisible(old)
}

set_male_light_moderate_barrier(10)
set_male_moderate_heavy_barrier(20)
set_male_binge_barrier(67.25)

#' Get and Set relative risk and prevalence/consumption data
#'


get_relative_risk_data <- function() {
  intermahpr_constants$relative_risk_data
}

set_relative_risk_data <- function(value) {
  old <- intermahpr_constants$relative_risk_data
  intermahpr_constants$relative_risk_data <- value
  invisible(old)
}

get_prevalence_consumption_data <- function() {
  intermahpr_constants$prevalence_consumption_data
}

set_prevalence_consumption_data <- function(value) {
  old <- intermahpr_constants$prevalence_consumption_data
  intermahpr_constants$prevalence_consumption_data <- value
  invisible(old)
}


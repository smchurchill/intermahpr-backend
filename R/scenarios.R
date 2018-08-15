#' Make a base interMAHP scenario
#' @export

makeNewModel <- function(rr, pc, mm) {
  free_rr <- rr %>%
    filterFree() %>%
    makeFreeFactories() %>%
    inner_join(pc, by = c("gender"))
  calibrated_rr <- rr %>%
    filterCalibrated() %>%
    makeCalibratedFactories(pc = pc, mm = mm)

  model <- bind_rows(free_rr, calibrated_rr) %>%
    select(getExpectedVars("model"))

  list(model = model, scenarios = list(), rr = rr, pc = pc, mm = mm)
  # %>% makeScenario(scenario_name = "Base", scale = 1)
}

#' make a scenario from a model object
#' @export

makeScenario <- function(.data, scenario_name = NA, scale) {
  pc <- rescale(.data = .data$pc, scale = scale) %>% computePopnMetrics()

  new_scenario <- .data$model %>%
    left_join(pc, by = c("region", "year", "gender", "age_group")) %>%
    mutate(
      fargs = pmap(
        list(
          p_fd = p_fd,
          mass = n_gamma,
          non_bingers = non_bingers,
          bingers = bingers,
          lb = lb,
          bb = bb,
          ub = ub),
        list
      )
    ) %>%
    mutate(
      current_fraction = map2(
        fargs,
        current_fraction_factory,
        ~.y(.x)
      )
    ) %>%
    mutate(
      former_fraction = map2(
        fargs,
        former_fraction_factory,
        ~.y(.x)
      )
    ) %>%
    select(getExpectedVars("scenario"))
  ## TODO: Add youngest x-y (into) 0-(x-1) agegroup.
  ## partially/former = 0,
  ## wholly = copy

  young <- queryYoung(.data)

  if(!is.null(young)) {
    young_scenario <- filter(new_scenario, age_group == young$young) %>%
      mutate(
        age_group = young$missing,
        current_fraction = map2(
          im, current_fraction,
          ~switch(
            checkAttributability(.x),
            "Partially" = function(...) 0,
            "Wholly" = .y
          )
        ),
        former_fraction = map(former_fraction, ~function(...) 0)
      )

    new_scenario <- bind_rows(young_scenario, new_scenario)
  }

  if(is.na(scenario_name)) scenario_name <- paste0("rescale_by_", scale)

  .data$scenarios[[scenario_name]] <- new_scenario

  .data
}

#' Get age-groups present in morbidity/mortality dataset but absent in prev-cons
#' dataset.
#'
#' @export
queryYoung <- function(.data) {
  mm_groups <- unique(.data$mm$age_group)
  pc_groups <- unique(.data$pc$age_group)

  missing <- setdiff(mm_groups, pc_groups)

  if(length(missing) == 0) {
    return(NULL)
  } else if(length(missing) > 1) {
    missing <- sort(missing)[1]
  }

  young <- sort(pc_groups)[1]

  list(missing = missing, young = young)
}

#' make multiple scenarios
#' @export

makeScenarios <- function(.data, scenario_names = NA, scales) {
  for(i in 1:length(scales)) {
    .data <- makeScenario(.data, scenario_names[i], scales[i])
  }
  .data
}

#' Compute a given scenario's AAF for former drinkers
#' @export
computeFormerFraction <- function(.data) {
  map_dbl(.data$former_fraction, ~.x())
}

#' Compute a given scenario's AAF for former drinkers and adds it to the
#' scenario
#' @export
addFormerFraction <- function(.data, var_name = "aaf_fd") {
  .data[[var_name]] <- computeFormerFraction(.data)
  .data
}

#' Compute a given scenario's AAF for current drinkers in the given intervals of
#' consumption, stratified over the given values of the gender variable, and add
#' it to the scenario.
#'
#'@param .data a scenario
#'@param strata a list where names(list) intersects with .data$gender, and each
#' entry of strata is a list with a "lower" and "upper" bound of consumption
#'@param var_name the name of the new variable to be added
#'@export
computeGenderStratifiedIntervalFraction <- function(.data, lower_strata, upper_strata) {
  .data %<>%
    mutate(
      x_lower = map_dbl(gender, ~`[[`(lower_strata, .x)),
      x_upper = map_dbl(gender, ~`[[`(upper_strata, .x)))

  .data$upper <- map2_dbl(.data$x_upper, .data$current_fraction , ~.y(.x))
  .data$lower <- map2_dbl(.data$x_lower, .data$current_fraction , ~.y(.x))
  .data$upper - .data$lower
}

#' Compute a given scenario's AAF for current drinkers in the given intervals of
#' consumption, stratified over the given values of the gender variable, and add
#' it to the scenario.
#'
#'@inheritParams computeGenderStratifiedIntervalFraction
#'@param var_name the name of the new variable to be added
#'@export
addGenderStratifiedIntervalFraction <- function(.data, lower_strata, upper_strata, var_name = "aaf_xd") {
  .data[[var_name]] <- computeGenderStratifiedIntervalFraction(.data, lower_strata, upper_strata)
}

#' Compute a given scenario's AAF for current drinkers in a given interval of
#' consumption
#' @export
computeIntervalFraction <- function(.data, lower = -Inf, upper = Inf) {
  map_dbl(.data$current_fraction, ~.x(upper) - .x(lower))
}

#' Compute a given scenario's AAF for current drinkers in a given interval of
#' consumption and adds it to the scenario
#' @export
addIntervalFraction <- function(.data, lower, upper, var_name = "aaf_xd") {
  .data[[var_name]] <- computeIntervalFraction(.data, lower, upper)
  .data
}

addIntervalFractions <- function(.data, lower, upper, grp_names) {
  n <- length(grp_names)
  if(!all(c(n == length(lower), n == length(upper)))) {
    warning("Groups not processed due to length mismatch")
    return(.data)
  }
  for(i in seq_len(n)) {
    .data <- addIntervalFraction(
      .data,
      lower[i],
      upper[i],
      paste0("aaf_", grp_names[i])
    )
  }
  .data
}

#' Compute a given scenario's AAF for current drinkers
#' @export
computeCurrentFraction <- function(.data) {
  computeIntervalFraction(.data)
}

#' Compute a given scenario's AAF for current drinkers and adds it to the
#' scenario
#' @export
addCurrentFraction <- function(.data, var_name = "aaf_cd") {
  .data[[var_name]] <- computeCurrentFraction(.data)
  .data
}

#' Compute a given scenario's Total AAF
#' @export
computeTotalFraction <- function(.data) {
  computeFormerFraction(.data) + computeCurrentFraction(.data)
}

#' Compute a given scenario's Total AAF and adds it to the scenario
#' @export
addTotalFraction <- function(.data, var_name = "aaf") {
  .data[[var_name]] <- computeTotalFraction(.data)
  .data
}

#' Extract and derive pertinent comparative scenario data from a model
#' @export
distillModel <- function(.data) {
  scenarios <- .data$scenarios
  master_name_list <- names(scenarios)
  aaf_name_list <- paste0("AAF: ", master_name_list)
  for(name in master_name_list) {
    scenario <- scenarios[[name]]
    scenario <- addTotalFraction(scenario, var_name = paste0("AAF: ", name)) %>%
      select(-contains("_fraction"))
    scenarios[[name]] <- scenario
  }

  by_vars <- getExpectedVars("distill_by")
  reduction <- reduce(scenarios, left_join, by = by_vars)

  attr <- (reduction$attributability == "Wholly")

  for(i in 2:length(master_name_list)) {
    reduction[[paste0("Relative AAF: ", master_name_list[i])]] <-
      reduction[[aaf_name_list[i]]] / reduction[["AAF: Base"]]

    reduction[[aaf_name_list[i]]] <-
      ifelse(attr, 1, reduction[[aaf_name_list[i]]])
  }

  reduction$`AAF: Base` <-
    ifelse(attr, 1, reduction$`AAF: Base`)

  reduction
}

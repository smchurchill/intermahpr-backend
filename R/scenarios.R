#' Make a base interMAHP scenario

makeNewModel <- function(rr, pc, dh) {
  free_rr <- rr %>%
    filterFree() %>%
    makeFreeFactories() %>%
    inner_join(pc, by = c("gender"))
  calibrated_rr <- rr %>%
    filterCalibrated() %>%
    makeCalibratedFactories(pc = pc, dh = dh)

  model <- bind_rows(free_rr, calibrated_rr) %>%
    select(getExpectedVars("model"))

  list(model = model, scenarios = list(), rr = rr, pc = pc, dh = dh) %>%
    makeScenario(scenario_name = "base", scale = 1)
}

#' make a scenario from a model object

makeScenario <- function(.data, scenario_name = NA, scale) {
  pc <- .data$pc
  if(scale != 1) {
    pc <- rescale(.data = pc, scale = scale) %>% computePopnMetrics()
  }

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

  if(is.na(scenario_name)) scenario_name <- paste0("rescale_by_", scale)

  .data$scenarios[[scenario_name]] <- new_scenario

  .data
}

#' make multiple scenarios

makeScenarios <- function(.data, scenario_names = NA, scales) {
  for(i in 1:length(scales)) {
    .data <- makeScenario(.data, scenario_names[i], scales[i])
  }
  .data
}

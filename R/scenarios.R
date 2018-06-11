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
    select(model_vars)

  list(model = model, scenarios = list(), rr = rr, pc = pc, dh = dh) %>%
    generateScenario(scenario_name = "base", scale = 1)
}

#' Generate a scenario from a model object

generateScenario <- function(.data, scenario_name = NA, scale) {
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
    select(scenario_vars)

  if(is.na(scenario_name)) scenario_name <- paste0("rescale_by_", scale)

  .data$scenarios[[scenario_name]] <- new_scenario

  .data
}

#' Generate multiple scenarios

generateScenarios <- function(.data, scenario_names = NA, scales) {
  for(i in 1:length(scales)) {
    .data <- generateScenario(.data, scenario_names[i], scales[i])
  }
  .data
}

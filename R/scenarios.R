#' Make a base interMAHP scenario

makeBaseScenario <- function(rr, pc, dh) {
  free_rr <- rr %>%
    filterFree() %>%
    makeFreeFactories() %>%
    inner_join(pc, by = c("gender"))
  calibrated_rr <- rr %>%
    filterCalibrated() %>%
    makeCalibratedFactories(pc = pc, dh = dh)

  bind_rows(free_rr, calibrated_rr) %>%
    select(
      incidence, region, year, gender, age_group, im, condition, outcome,
      current_fraction_factory, former_fraction_factory)
}

makeScaledScenario <- function(scenario, scale) {



}




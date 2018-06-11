#### Population Specific Data Carpentry ----------------------------------------

#' Prepare Population Data

preparePC <- function(.data, ...) {
  pc %<>%
    cleanPC() %>%
    setPopnConstants(...) %>%
    computePopnMetrics()
}

#' Clean Population data

cleanPC <- function(.data) {
  clean(pc, pc_vars)
}

#' List of variables expected to be in an input PC sheet
#'
#'

pc_vars <- c(
  "region",
  "year",
  "gender",
  "age_group",
  "population",
  "pcc_litres_year",
  "correction_factor",
  "relative_consumption",
  "p_la",
  "p_fd",
  "p_cd",
  "p_bd"
)


#' List of pc constants
#'
#'

pc_constants <- c(
  "bb",
  "lb",
  "ub"
)

#### Population metric alterations ---------------------------------------------

#' Set Population Constants
#'@param bb binge barrier
#'@param lb lower bound of consumption
#'@param ub upper bound of consumption
#'

setPopnConstants <- function(
  .data, bb = list("Female" = 53.8, "Male" = 67.25), lb = 0.03, ub = 250
) {
  mutate(.data, lb = lb, bb = map_dbl(gender, ~`[[`(bb, .x)), ub = ub)
}

#' Compute Population Metrics

computePopnMetrics <- function(.data) {
  ## Magic numbers
  gc = list("Female" = 1.582564, "Male" = 1.371241)
  yearly_to_daily_conv = 0.002739726
  litres_to_millilitres_conv = 1000
  millilitres_to_grams_ethanol_conv = 0.7893

  .data %>%
    group_by(region, year) %>%
    mutate(
      pcc_g_day =
        pcc_litres_year *
        litres_to_millilitres_conv *
        millilitres_to_grams_ethanol_conv *
        yearly_to_daily_conv *
        correction_factor,
      drinkers = population * p_cd,
      pcad = pcc_g_day * sum(population) / sum(drinkers),
      pcc_among_drinkers = relative_consumption * pcad * sum(drinkers) /
        sum(relative_consumption*drinkers)
    )%>%
    ungroup %>%
    mutate(
      gamma_constant = map_dbl(gender, ~`[[`(gc, .x))
    ) %>%
    mutate(
      gamma_shape = 1/gamma_constant,
      gamma_scale = gamma_constant*pcc_among_drinkers
    ) %>%
    mutate(
      glb = pgamma(q = lb, shape = gamma_shape, scale = gamma_scale),
      gbb = pgamma(q = bb, shape = gamma_shape, scale = gamma_scale),
      gub = pgamma(q = ub, shape = gamma_shape, scale = gamma_scale)
    ) %>%
    mutate(
      nc = gub - glb
    ) %>%
    mutate(
      df = p_cd / nc
    ) %>%
    mutate(
      n_gamma = pmap(
        list(shape = gamma_shape, scale = gamma_scale, factor = df),
        makeNormalizedGamma
      ),
      p_bat = df * (gub - gbb)
    ) %>%
    mutate(
      n_pgamma = pmap(list(f = n_gamma, lb = lb), makeIntegrator)
    ) %>%
    mutate(
      non_bingers = (p_cd - p_bd)  / (p_cd - p_bat),
      bingers = (p_bd - p_bat) / (p_cd - p_bat)
    )
}

#' Rescale Population Data
#'
#'@param .data cleaned population data with constants set
#'@param scale a percentage of the current consumption expected in the
#'scenario under study
#'
#'@return Rescaled consumption data (just pc_vars)
#'

rescale <- function(.data, scale = 1) {
  base <- computePopnMetrics(.data)

  .data %>%
    mutate(pcc_litres_year = scale * pcc_litres_year) %>%
    computePopnMetrics() %>%
    mutate(p_bd = p_bd * p_bat / base$p_bat) %>%
    select(c(pc_vars, pc_constants))
}

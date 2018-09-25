## intermahpr - R package backend for the intermahp shiny app
## Copyright (C) 2018 Canadian Institute for Substance Use Research

#### Population Specific Data Carpentry ----------------------------------------

#' Prepare Population Data
#' @export
preparePC <- function(.data, ...) {
  message("Preparing prevalence and consumption input... ", appendLF = FALSE)
  .data %<>%
    clean(getExpectedVars("pc")) %>%
    setPopnConstants(...) %>%
    computePopnMetrics()

  message("Done")
  .data
}

#' Return PC dataset for viewing in wide format
#' @export
renderPCWide <- function(.data) {
  .data %>%
    select(getExpectedVars("pc_display")) %>%
    rename("gamma_normalizer" = "nc")
}

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
  ## 'Magic' numbers
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
      ## alcohol consumption over all age groups
      pcad = pcc_g_day * sum(population) / sum(drinkers),
      ## mean consumption per age group
      pcc_among_drinkers = relative_consumption * pcad * sum(drinkers) /
        sum(relative_consumption*drinkers)
    )%>%
    ungroup %>%
    mutate(
      gamma_shape = 1/gamma_constant/gamma_constant,
      gamma_scale = gamma_constant*gamma_constant*pcc_among_drinkers
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
    select(getExpectedVars("pc", "constants"))
}

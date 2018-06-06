compute_pc <- function(
  .data,
  bb = list("Female" = 53.8, "Male" = 67.25),
  lb = 0.03,
  ub = 250
) {
  ## Magic numbers
  gc = list("Female" = 1.582564, "Male" = 1.371241)
  yearly_to_daily_conv = 0.002739726
  litres_to_millilitres_conv = 1000
  millilitres_to_grams_ethanol_conv = 0.7893

  .data %>%
    group_by(REGION, YEAR) %>%
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
      gamma_scale = gamma_constant*pcc_among_drinkers,
      bb = map_dbl(gender, ~`[[`(bb, .x)),
      lb = lb,
      ub = ub
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
        normalized_gamma_factory
      ),
      p_bat = df * (gub - gbb)
    ) %>%
    mutate(
      r1 = (p_cd - p_bd)  / (p_cd - p_bat),
      r2 = (p_bd - p_bat) / (p_cd - p_bat)
    )
}



#' Scales the per capita consumption and binge drinker prevalence
#'
#'@param scale a percentage of the current consumption expected in the
#'scenario under study
#'
#'

scale_pc <- function(
  .data,
  scale = 1,
  bb = list("Female" = 53.8, "Male" = 67.25),
  lb = 0.03,
  ub = 250
) {



  base_f <- format_pc(.data)
  name_r <- names(base_f)
  base_d <- derive_pc(base_f)

  .data %>%
    mutate(PCC_litres_year = scale * PCC_litres_year) %>%
    format_pc() %>%
    derive_pc() %>%
    mutate(P_BD = P_BD * P_BAT / base_d$P_BAT) %>%
    select(name_r)
}


#' List of variables expected to be in a PC sheet
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

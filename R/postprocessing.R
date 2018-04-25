#' Postprocessing function that tidies base AAF into a certain outcome type
#'
#'@description
#'  Filters AAFs into any of mortality, morbidity, or combined.  Exported for
#'  use in the Shiny app, called behind the scenes when invoking intermahp_base.
#'
#'@importFrom magrittr "%>%" "%<>%"
#'@importFrom dplyr filter
#'
#'@export
#'


outcome_splitter <- function(aaf_table, outcome) {
  aaf_table %<>%
    filter(OUTCOME == outcome | OUTCOME == "Combined") %>%
    mutate(AAF_LD = sapply(AAF_GRP, `[[`, 1),
           AAF_MD = sapply(AAF_GRP, `[[`, 2),
           AAF_HD = sapply(AAF_GRP, `[[`, 3)) %>%
    select(REGION, YEAR, GENDER, AGE_GROUP, IM, CONDITION,
           AAF_FD, AAF_LD, AAF_MD, AAF_HD, AAF_TOTAL)

  aaf_table
}


#' Postprocessing function that integrates a list of functions from a list of
#' lower bounds to a list of upper bounds.
#'
#' Integrate that takes as input a list funs of functions, a vector vlower of
#' "lower" values and a vector vupper of "upper" values, assumes that
#' length(vlower) = length(vupper) = length(funs), and returns a
#' vector of integrals of funs between the lower and upper values.
#'
#'
#'

vintegrate <- function(funs, vlower, vupper) {
  llower <- length(vlower)
  lupper <- length(vupper)
  lfuns  <- length(funs)

  values <- rep(-1, lfuns)

  for(i in 1:lfuns) {
    values[[i]] <- integrate(funs[[i]], vlower[[i]], vupper[[i]])$value
  }

  values
}

#' Helper function that tidies base AAF into the relevant prevalence and
#' consumption data
#'
#'@description
#' Intended for use with intermahp base functionality (assumes the positions of
#' cutpoint variables).  Exported for use in the Shiny app, called behind the
#' scenes when invoking intermahp_base.
#'
#'@param aaf_table a tibble as returned by add_cutpoints with cutpoints of the
#'form c(LB, LM, MH, UB) where LM is the light-moderate barrier and MH is the
#'moderate-heavy barrier.
#'
#'@importFrom dplyr distinct select case_when
#'
#'@export
#'

extract_prevcons <- function(aaf_table) {
  aaf_table %<>%
    mutate(
      LM = sapply(CUTS, `[[`, 2),
      MH = sapply(CUTS, `[[`, 3)
    ) %>%
    distinct(
      REGION, YEAR, GENDER, AGE_GROUP, POPULATION, PCC_AMONG_DRINKERS,
      GAMMA_SHAPE, GAMMA_SCALE, P_LA, P_FD, P_CD, EXT, LB, LM, MH, UB,
      .keep_all = TRUE
    ) %>%
    select(
      REGION, YEAR, GENDER, AGE_GROUP, POPULATION, PCC_AMONG_DRINKERS,
      GAMMA_SHAPE, GAMMA_SCALE, P_LA, P_FD, P_CD, EXT, LB, LM, MH, UB, N_GAMMA
    )

  aaf_table %<>%
    mutate(
      EXTRAPOLATION = case_when(
        EXT == TRUE  ~ "linear",
        EXT == FALSE ~ "capped")
    ) %>%
    add_column(
      P_LD = vintegrate(
        aaf_table[["N_GAMMA"]],
        aaf_table[["LB"]],
        aaf_table[["LM"]]),
      P_MD = vintegrate(
        aaf_table[["N_GAMMA"]],
        aaf_table[["LM"]],
        aaf_table[["MH"]]),
      P_HD = vintegrate(
        aaf_table[["N_GAMMA"]],
        aaf_table[["MH"]],
        aaf_table[["UB"]])
    ) %>%
    mutate(P_CD_SUM = P_LD + P_MD + P_HD) %>%
    select(
      REGION, YEAR, GENDER, AGE_GROUP, POPULATION, PCC_AMONG_DRINKERS,
      GAMMA_SHAPE, GAMMA_SCALE, P_LA, P_FD, P_CD, P_LD, P_MD, P_HD,
      EXTRAPOLATION
    )

  aaf_table
}

#' Compute Alcohol Attributable Counts
#'
#'@description
#'  Applies the computed AAFs to user provided death/hospitalization counts.
#'
#'@param aaf_table a tibble of AAFs as returned by compute_aafs.
#'@param dh a formatted deaths/hosps table as returned by format_v*_dh
#'
#'

aa_counts <- function(aaf_table, dh) {
}

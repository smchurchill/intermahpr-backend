##### g02-postprocessors #######################################################
##
## Postprocessing functions called from functions in s02-process
##
##
##
##

#' Compute AAFs for all conditions for current drinkers and totals AAFs
#'
#'@param aaf_table A tibble as returned by assemble (i.e. has the following vars
#' AAF_CMP: fn, VIVO, computes AAF from LB to each input
#' AAF_FD: dbl, AAF for former drinkers
#' UB: dbl, consumption upper bound
#')
#'
#'@importFrom purrr pmap map map2
#'@importFrom magrittr "%>%" "%<>%"
#'
#'@export
#'

aaf_total <- function(aaf_table) {
  aaf_table %>%
    mutate(AAF_CD = map2_dbl(AAF_CMP, UB, ~(.x(.y)))) %>%
    mutate(AAF_TOTAL = AAF_FD + AAF_CD)
}


#' Compute AAFs for all conditions as separated by the given cutpoints
#'
#'@param aaf_table A tibble as returned by assemble (i.e. has the following vars
#' GENDER: chr, gender levels
#' AAF_CMP: fn, VIVO, computes AAF from LB to each input
#' AAF_FD: dbl, AAF for former drinkers
#' LB: dbl, consumption lower bound
#' UB: dbl, consumption upper bound
#')
#'@param cuts A sorted list of double vectors indexed by aaf_table$GENDER st.
#' each value is between LB and UB
#'
#'@importFrom purrr pmap map map2
#'@importFrom magrittr "%>%" "%<>%"
#'
#'@export
#'

compute_aafs <- function(aaf_table, cuts) {
  aaf_table %>%
    mutate(CUTS = pmap(list(LB, cuts[GENDER], UB), c)) %>%
    mutate(CUMUL_F = map2(AAF_CMP, CUTS, ~.x(.y))) %>%
    mutate(AAF_GRP = map(CUMUL_F, ~diff(.x))) %>%
    mutate(AAF_CD = map2_dbl(AAF_CMP, UB, ~(.x(.y)))) %>%
    mutate(AAF_TOTAL = AAF_FD + AAF_CD)
}

#' Add evaluation cutpoints to a given datatable
#'
#' Written to provide additional reactivity to Shiny InterMAHP app.
#'
#'@param aaf_table a tibble as returned by assemble
#'@param cuts a list of double vectors indexed by aaf_table$GENDER
#'
#'@importFrom purrr pmap
#'@importFrom magrittr "%<>%"
#'
#'@export
#'

add_cutpoints <- function(aaf_table, cuts) {
  aaf_table %>%
    mutate(CUTS = pmap(list(LB, cuts[GENDER], UB), c))
}

#' Evaluate a given datatable at pre-added cutpoints
#'
#' Written to provide additional reactivity to shiny app
#'
#'@param aaf_table_cuts a tibble as returned by add_cutpoints
#'
#'@importFrom purrr pmap map map2
#'@importFrom magrittr "%>%" "%<>%"
#'
#'@export
#'

evaluate_at_cutpoints <- function(aaf_table_cuts) {
  aaf_table_cuts %>%
    mutate(CUMUL_F = map2(AAF_CMP, CUTS, ~.x(.y))) %>%
    mutate(AAF_GRP = map(CUMUL_F, ~diff(.x))) %>%
    mutate(AAF_CD = map2_dbl(AAF_CMP, UB, ~(.x(.y)))) %>%
    mutate(AAF_TOTAL = AAF_FD + AAF_CD)
}


#' Postprocessing function that tidies base AAF by outcome type
#'
#'@description
#'  Splits "Combined" levels of an OUTCOME variable into Morbidity and Mortality
#'
#'@param .data Any table with an OUTCOME variable whose levels are:
#'  "Mortality", "Morbidity", and "Combined".
#'
#'@return .data with an OUTCOME variable whose levels are: "Mortality" and
#'  "Morbidity" and 2c+t+b rows, where c is the number of OUTCOME = "Combined"
#'  rows in aaf_table, t is the number of OUTCOME = "Mortaility", and b is the
#'  number of OUTCOME  = "Morbidity".
#'
#'@importFrom magrittr "%>%" "%<>%"
#'@importFrom dplyr filter
#'
#'@export
#'

split_outcome <- function(.data) {
  mortality <- filter(.data, OUTCOME == "Mortality" | OUTCOME == "Combined")
  morbidity <- filter(.data, OUTCOME == "Morbidity" | OUTCOME == "Combined")

  mortality$OUTCOME <- "Mortality"
  morbidity$OUTCOME <- "Morbidity"

  rbind(mortality, morbidity)
}



#' Postprocessing function that extracts Light, Moderate, Heavy group AAFs
#'
#'
#'
#'@export
#'

name_cuts <- function(aaf_table) {
  aaf_table %<>%
    mutate(
      AAF_LD = vapply(AAF_GRP, `[[`, 1, FUN.VALUE = 0),
      AAF_MD = vapply(AAF_GRP, `[[`, 2, FUN.VALUE = 0),
      AAF_HD = vapply(AAF_GRP, `[[`, 3, FUN.VALUE = 0)) %>%
    select(
      REGION, YEAR, GENDER, AGE_GROUP, IM, CONDITION, OUTCOME,
      AAF_FD, AAF_LD, AAF_MD, AAF_HD, AAF_TOTAL)
  aaf_table
}


#' Postprocessing helper function that integrates a list of functions from a
#' list of lower bounds to a list of upper bounds.
#'
#' Integrate that takes as input a list funs of functions, a vector vlower of
#' "lower" values and a vector vupper of "upper" values, assumes that
#' length(vlower) = length(vupper) = length(funs), and returns a
#' vector of integrals of funs between the lower and upper values.
#'
#' This function is not optimized or "vectorized," it simply assigns memory and
#' uses a for loop to perform each integral.
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
      LM = vapply(CUTS, `[[`, 2, FUN.VALUE = 0),
      MH = vapply(CUTS, `[[`, 3, FUN.VALUE = 0)
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

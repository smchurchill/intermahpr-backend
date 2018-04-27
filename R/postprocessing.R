#' Postprocessing function that tidies base AAF by outcome type
#'
#'@description
#'  Splits "Combined" levels of an OUTCOME variable into Morbidity and Mortality
#'
#'@param aaf_table Any table with an OUTCOME variable whose levels are:
#'  "Mortality", "Morbidity", and "Combined".
#'
#'@return aaf_table with an OUTCOME variable whose levels are: "Mortality" and
#'  "Morbidity" and 2c+t+b rows, where c is the number of OUTCOME = "Combined"
#'  rows in aaf_table, t is the number of OUTCOME = "Mortaility", and b is the
#'  number of OUTCOME  = "Morbidity".
#'
#'@importFrom magrittr "%>%" "%<>%"
#'@importFrom dplyr filter
#'
#'@export
#'

outcome_splitter <- function(aaf_table) {
  mortality <- filter(aaf_table, OUTCOME == "Mortality" | OUTCOME == "Combined")
  morbidity <- filter(aaf_table, OUTCOME == "Morbidity" | OUTCOME == "Combined")

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
      AAF_LD = sapply(AAF_GRP, `[[`, 1),
      AAF_MD = sapply(AAF_GRP, `[[`, 2),
      AAF_HD = sapply(AAF_GRP, `[[`, 3)) %>%
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


## Split Similar AAF rows ----
#' Split Similar AAF rows
#'
#'@description
#'  condition specific method.
#'   - 4.5 AAFs are proportionally distributed to 4.4, 4.6, and 4.7.
#'   - 5.3 receives the AAF distribution of LD: 0, MD: 0, HD: 0.
#'   - 8.all and 9.all are split exactly over 8.1-4, 8.6 and 9.1, 9.3-5 resp.,
#'       and proportionally over 8.5 and 9.2 resp.
#'
#'@param aaf_table as produced by name_cuts (i.e has names(aaf_table) of
#'  REGION, YEAR, GENDER, AGE_GROUP, IM, CONDITION, OUTCOME,
#'  AAF_FD, AAF_LD, AAF_MD, AAF_HD, AAF_TOTAL)
#'
#'

split_sim <- function(aaf_table) {
  epi_sim <- epilepsy <- filter(aaf_table, IM == "(4).(5)")
  epi_aaf <- epilepsy[c("AAF_LD", "AAF_MD", "AAF_HD", "AAF_TOTAL")]

  epi_xaaf <- epi_aaf / epi_aaf$AAF_TOTAL
  epi_sim[c("AAF_LD", "AAF_MD", "AAF_HD", "AAF_TOTAL")] <- epi_xaaf

  c44 <- "Degeneration of nervous system due to alcohol"
  c46 <- "Alcoholic polyneuropathy"
  c47 <- "Alcoholic myopathy"

  degen <- mutate(epi_sim, IM = "(4).(4)", CONDITION = c44)
  polyn <- mutate(epi_sim, IM = "(4).(6)", CONDITION = c46)
  myopa <- mutate(epi_sim, IM = "(4).(7)", CONDITION = c47)






}





## Compute Alcohol Attributable Counts ----
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

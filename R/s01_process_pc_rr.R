#### Process Derived Data ------------------------------------------------------

#' Joins PC and RR and obtains AAF functions
#'
#'@description Performs a full_join of RR and PC over GENDER and then obtains an
#'  AAF function and AAF_FD computation for each observation.  Does a sanity
#'  check, ensuring that the joining variable is well-matched.
#'
#'@param rr relative risk tibble as produced by derive_v*_rr
#'@param pc prevalence and consumption tibble as produced by derive_v*_pc
#'
#'@return tibble with one row per unique region.year.gender.age_group.im combn,
#'  an AAF_FD variable with alc.attr. fraction for former drinkers, and an
#'  AAF_CMP fn and aaf_cd, aaf_total dbl variables
#'
#'@export
#'


join_pc_rr <- function(pc, rr) {
  JOINT <- dplyr::full_join(pc, rr, by = "GENDER")

  if(any(is.na(JOINT))) {
    MANGLED <- JOINT[rowSums(is.na(JOINT)) > 0,]
    IGNORE <- MANGLED[, c("REGION", "YEAR", "GENDER",
                        "AGE_GROUP", "IM", "CONDITION")]

    ignore_message <- paste0(capture.output(IGNORE), collapse = "\n")

    warning(
      paste0(
        "\nNAs introduced by joining of prevalence and consumption data with ",
        "relative risk data.\nTo avoid this, ensure that the levels in the ",
        "Gender column of prevalence and consumption input match the levels ",
        "in the Gender column of the relative risk input exactly.\nThe ",
        "following observations will be ignored:\n",
        collapse = "\n"
      ),
      ignore_message
    )
    JOINT <- JOINT[rowSums(is.na(JOINT)) == 0,]
  }

  JOINT %>%
    mutate(
      INTGRND = pmap(
        list(BB, LB, UB, R1, R2, LNXT_RR, BNGD_RR, N_GAMMA),
        intgrnd_factory
      )
    ) %>%
    mutate(
      AAF_CMP = pmap(
        list(LB, UB, RR_FD, P_FD, INTGRND),
        aaf_cmp_factory
      )
    ) %>%
    mutate(
      AAF_FD = pmap_dbl(
        list(LB, UB, RR_FD, P_FD, INTGRND),
        aaf_fd
      ),
      AAF_CD = map2_dbl(AAF_CMP, UB, ~(.x(.y)))
    ) %>%
    mutate(
      AAF_TOTAL = AAF_CD + AAF_FD
    )
}

#' Collect and assemble AAF data from formatted RR and PC data
#'
#'@param pc  Prevalence / Consumption input as produced by format_v*_pc
#'@param rr  Relative Risk input as produced by format_v*_rr
#'@param ext logical, extrapolate linearly?
#'@param lb  Double, consumption lower bound
#'@param ub  Double, consumption upper bound
#'@param bb  Double vector, Binge consumption level, Gender stratified
#'@param gc  Gamma constant.  The linear relationship between mean and standard
#'  deviation within the gamma distribution that describes consumption among
#'  current drinkers.  Stratified by gender -- names(gc) must match levels of
#'  rr$GENDER and pc$GENDER.
#'
#'

assemble <- function(pc, rr, ext, lb, ub, bb, gc) {
  RRD <- derive_v0_rr(rr = rr, ext = ext)
  PCD <- derive_v0_pc(pc = pc, bb = bb, lb = lb, ub = ub, gc = gc)

  join_pc_rr(pc = PCD, rr = RRD)
}




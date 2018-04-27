#### Process Derived Data ------------------------------------------------------

#' Joins PC and RR and obtains AAF functions
#'
#'@description Performs a full_join of RR and PC over GENDER and then obtains an
#'  AAF function and AAF_FD computation for each observation.  Does a sanity
#'  check, ensuring that the joining variable is well-matched.
#'
#'@param RR relative risk tibble as produced by derive_v*_rr
#'@param PC prevalence and consumption tibble as produced by derive_v*_pc
#'
#'@return tibble with one row per unique region.year.gender.age_group.im combn,
#'  an AAF_FD variable with alc.attr. fraction for former drinkers, and an
#'  AAF_CMP variable that
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

  for(n in 1:nrow(JOINT)) {
    intgrnd <- intgrnd_factory(JOINT[n, ])
    JOINT[[n, "INTGRND"]] <- intgrnd

    aaf_cmp <- aaf_cmp_factory(JOINT[n, ])
    JOINT[[n, "AAF_CMP"]] <- aaf_cmp

    JOINT[[n, "AAF_FD"]] <- aaf_fd(JOINT[n, ])
  }

  JOINT
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

#' Combine Prev/Cons and Death/Hosp data to calibrate wholly attributable curves
#'
#'@description
#' Wholly attributable conditions have an AAF_TOTAL of 1.00, but we can still
#' distribute this mass over the interval of interest (i.e. 0.03 to UB).
#' When applicable (4.(1, 2, 3), 6.(1, 5) we calibrate a loglinear conditional
#' probability curve whose area is equal to the regional yearly prevalence of
#' the given condition among current drinkers.
#'
#'@param dh as returned by derive_v*_dh
#'
#'@importFrom magrittr %<>% %>%
#'

calibrate <- function(dh) {
  ## Applicable IMs
  CABLE <- c("(4).(1)", "(4).(2)", "(4).(3)", "(6).(1)", "(6).(5)")

  dh %<>%
    filter(IM %in% CABLE)


  ## Start with a tibble of wholly attributable conditions and how we deal with
  ## them
  WA <- tibble::tibble(
    IM = c(
      "(4).(1)",
      "(4).(2)",
      "(4).(3)",
      "(6).(1)",
      "(6).(5)"
    ),
    THRESHOLD = c(
      "BB",
      "BB",
      "BB",
      "LB",
      "LB"
    )
  )
}

#' Classify IMs for AAF treatment
#'
#'
#'

classify <- function(dh) {
  ## Calibrate
  CABLE <- tibble::tibble(
    IM =    list("(4).(1)", "(4).(2)", "(4).(3)", "(6).(1)", "(6).(5)"),
    USING = list("BB",      "BB",      "BB",      "LB",      "LB")
  )

  ## Distribute proportionally
  DISTR <- tibble::tibble(
    IM =    list("(4).(4)", "(4).(6)", "(4).(7)", "(8).(5)",   "(9).(2)"),
    USING = list("(4).(5)", "(4).(5)", "(4).(5)", "(8).(all)", "(9).(all)")
    )

  ## Copy directly
  COPY <- tibble::tibble(
    IM =    list(
      "(5).(7)",
      "(8).(1)",   "(8).(2)",   "(8).(3)",   "(8).(4)",   "(8).(6)",
      "(9).(1)",   "(9).(3)",   "(9).(4)",   "(9).(5)"
    ),
    USING = list(
      "(6).(2)",
      "(8).(all)", "(8).(all)", "(8).(all)", "(8).(all)", "(8).(all)",
      "(9).(all)", "(9).(all)", "(9).(all)", "(9).(all)"
    )
  )
}




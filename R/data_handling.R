#### Handle missing data -------------------------------------------------------

#' What value is imputed when a variable is missing?
#'
#' In order to normalize code in future methods, missing data is removed either
#' through imputation when possible, or by throwing an error when the data is
#' necessary.  impute_with is called when a missing value \emph{can} be imputed.
#'
#'
#' @param var is a string -- a variable name.
#'
#' @return A string or integer that is appropriate for imputation.
#'
#'

impute_with <- function(var) {
  if(grepl("B[0-9]", var)) {var = "BETA"}
  switch(
    var,
    IM = "(0).(0)",
    CONDITION = "Unspecified",
    REGION = "Unspecified",
    YEAR = "Unspecified",
    AGE_GROUP = "Unspecified",
    OUTCOME = "Combined",
    FUNCTION = "FP",
    RR_FD = 1,
    BINGEF = 1,
    BETA = 0,
    COUNT = 0
  )
}

#' Definition of missing data for the purposes of intermahpr
#'

is.missing <- function(obs) {
  is.na(obs) | is.null(obs) | obs == "."
}

#' Dummy function for allocating memory
#'
#'@param ... Accepts any input
#'
#'@return Always returns 0

zero <- list(fn = function(...) 0)


#### Verify and wrangle data ---------------------------------------------------

#' Format relative risk input data to our desired specifications
#'
#'@description
#'Formats all variable names to upper case and imputes imputable missing data.
#'Intended for a relative risk sheet as described in the intermahp user guide.
#'
#'expected variables:
#'  Gender (paired with PC, DH data)
#'  Function (FP, Spline, with valid IM, or Step, for HIV only)
#'  B1-B16 (only if Function == FP. need nonzero betas, rest can be missing)
#'  IM
#'  Condition
#'  Outcome
#'  RR_FD
#'  BingeF
#'
#'@param rr is a dataset that can be converted to tibble with the expected
#'  variables
#'
#'@export
#'

format_v0_rr <- function(rr) {
  RR <- tibble::as.tibble(rr)
  names(RR) <- do.call(stringr::str_to_upper, list(names(RR)))
  BETAS <- do.call(paste0, list(rep("B", 16), 1:16))
  EXPECTED <- c(
    "IM", "CONDITION", "GENDER", "OUTCOME",
    "RR_FD", "BINGEF", "FUNCTION", BETAS
  )
  MISSING <- EXPECTED[!(EXPECTED %in% names(RR))]
  if(length(MISSING) > 0) {
    stop(
      "Missing variables from the supplied relative risk sheet: ",
      paste(MISSING, collapse = ", ")
    )
  }

  for(var in names(RR)) {
    RR[var][is.missing(RR[var])] <- impute_with(var)
  }

  readr::type_convert(RR[EXPECTED])
}

#' Format prevalence and consumption input data to our desired specifications
#'
#'@description
#'Formats all variable names to upper case and imputes imputable missing data.
#'Intended for a prevalence and consumption sheet as described in the intermahp
#'user guide.
#'
#' expected variables:
#'    P_FD (Proportion, Former Drinkers)
#'    P_BD (Proportion, Binge Drinkers)
#'    P_LA (Proportion, Lifetime Abstainers)
#'    P_CD (Proportion, Current Drinkers)
#'    Gamma Distribution Parameters:
#'      PCC_litres_year (Per Capita Consumption)
#'      Population
#'      Relative_Consumption
#'      Correction_Factor
#'    Gender (in order to pair with RR curves)
#'    Region
#'    Year
#'    Gender
#'    Age_Group
#'
#'@param rr is a dataset that can be converted to tibble with the expected
#'  variables
#'
#'@export
#'

format_v0_pc <- function(pc) {
  PC <- tibble::as.tibble(pc)
  names(PC) <- do.call(stringr::str_to_upper, list(names(PC)))


  EXPECTED <- c(
    "REGION", "YEAR", "GENDER", "AGE_GROUP", "POPULATION", "PCC_LITRES_YEAR",
    "CORRECTION_FACTOR", "RELATIVE_CONSUMPTION", "P_LA", "P_FD", "P_CD", "P_BD"
  )
  MISSING <- EXPECTED[!(EXPECTED %in% names(PC))]
  if(length(MISSING) > 0) {
    stop(
      "Missing variables from the supplied prevalence/consumption sheet: ",
      paste(MISSING, collapse = ", ")
    )
  }

  for(var in names(PC)) {
    PC[var][is.missing(PC[var])] <- impute_with(var)
  }

  readr::type_convert(PC[EXPECTED])
}

#' Format Count Data
#'
#'@description
#'Formats all variable names to upper case and imputes imputable missing data.
#'Intended for a count sheet
#'
#' expected variables:
#'    IM (Used to join this data with aafs)
#'    Region
#'    Year
#'    Gender
#'    Age_Group
#'    Outcome
#'    Count
#'
#'@param dh is a dataset that can be converted to tibble with the expected
#'  variables
#'
#'@export
#'

format_v0_dh <- function(dh) {
  DH <- tibble::as.tibble(dh)
  names(DH) <- do.call(stringr::str_to_upper, list(names(DH)))

  EXPECTED <- c(
    "IM", "CONDITION", "REGION", "YEAR", "GENDER", "AGE_GROUP",
    "OUTCOME", "COUNT"
  )
  MISSING <- EXPECTED[!(EXPECTED %in% names(DH))]
  if(length(MISSING) > 0) {
    stop(
      paste0(
        "Missing and needed variables from the supplied counts sheet: "
      ),
      paste(MISSING, collapse = ", ")
    )
  }

  for(var in names(DH)) {
    DH[var][is.missing(DH[var])] <- impute_with(var)
  }

  readr::type_convert(DH[EXPECTED])
}

#### Process formatted data ----------------------------------------------------

#' Derives relative risk curves from relative risk data
#'
#'@param rr tibble as formatted above
#'@param ext logical indicator of extrapolation, (TRUE = linear, FALSE = capped)
#'
#'@return tibble which is RR with an additional column that contains all
#'  relevent relative risk curves under the variables BASE_RR, LNXT_RR, and
#'  BNGD_RR.
#'    RR[["BASE_RR"]] is the base relative risk function (i.e. no extrapolation
#'      after 100/150)
#'    RR[["LNXT_RR"]] is the extrapolated relative risk curve
#'    RR[["BNGD_RR"]] is the extrapolated relative risk curve for binge drinkers
#'
#'@export
#'

derive_v0_rr <- function(rr, ext) {
  rr[, "EXT"] <- ext
  rr <- add_column(
    rr,
    BASE_RR = zero,
    LNXT_RR = zero,
    BNGD_RR = zero
  )

  for(n in 1:nrow(rr)) {
    base <- set_rr(rr[n,])
    rr[[n, "BASE_RR"]] <- base

    lnxt <- ext_rr(rr[n,])
    rr[[n, "LNXT_RR"]] <- lnxt

    bngd <- bng_rr(rr[n,])
    rr[[n, "BNGD_RR"]] <- bngd
  }

  rr
}


#' Derives probability distribution parameters from Prevalence and Consumption
#' data
#'
#' From a given prevalence and consumption dataframe (or similar), computes re-
#' -quired data for AAF computation including Per Capita Consumption among
#' drinkers, and Gamma Shape/Scale parameters.  Also computes constants
#' necessary for AAF integration, namely the Deflation factor and Binge-Split
#' ratios.
#'
#' Gamma parameters are derivable entirely from prevalence and consumption,
#' no additional input is necessary.
#'
#'@param pc is assumed type tibble with names(PC) = c(YEAR, REGION, GENDER
#' , AGE_GROUP, POPULATION, PCC_LITRES_YEAR, CORRECTION_FACTOR,
#' RELATIVE_CONSUMPTION, P_LA, P_FD, P_CD, P_BD).
#'
#' Note that the bundled data set pc_default satisfies these constraints.
#'
#'@param bb User supplied gender stratified binge barrier
#'@param lb Lower bound
#'@param ub Upper bound
#'@param gc Gender stratified gamma constants
#'
#' Constants used:
#'  0.002739726 = 1/365
#'  0.7893 = unit conversion (millilitres to grams) for ethanol
#'
#'@return A tibble with additional variables:
#'   PCC_G_DAY: Per capita consumption, converted from litres/year to grams/day
#'   DRINKERS: Total drinkers in region, computed as POPULATION*P_CD
#'   PCAD: Average consumption among drinkers
#'   PCC_AMONG_DRINKERS: Computed from PCAD and relative consumption
#'   GAMMA_CONSTANT: The distribution of drinkers at drinking level of x grams
#'     per day is fit by a gamma distribution whose mean and standard deviation
#'     are related by GAMMA_CONSTANT, stratified by gender
#'   GAMMA_SHAPE: Computed from GAMMA_CONSTANT
#'   GAMMA_SCALE: Computed from GAMMA_CONSTANT
#'   INDEX: Among vector inputs, which index refers to this cohort?
#'   BB: Binge barrier
#'   LB: Lower bound of integration (0.03)
#'   UB: Upper bound of integration (user-specified, typical values are 150,250)
#'   NC: Normalizing constant for the gamma distribution (how much of the gamma
#'   distribution is between LB and UB?)
#'   DF: Deflation factor for the gamma distribution (DF = P_CD/NC)
#'   P_BAT: Proportion of bingers above threshold (BB)
#'   R1, R2: Ratios computed at this step and used in AAF computations
#'
#'@importFrom magrittr "%>%" "%<>%"
#'@importFrom dplyr group_by mutate
#'@importFrom tibble add_column
#'
#'@export
#'

derive_v0_pc <- function(pc, bb, lb, ub, gc) {
  PC <- pc
  PC %<>%
    group_by(REGION, YEAR) %>%
    mutate(
      PCC_G_DAY = PCC_LITRES_YEAR * 1000 *
        0.7893 * CORRECTION_FACTOR * 0.002739726,
      DRINKERS = POPULATION * P_CD,
      PCAD = PCC_G_DAY * sum(POPULATION) / sum(DRINKERS),
      PCC_AMONG_DRINKERS = RELATIVE_CONSUMPTION * PCAD * sum(DRINKERS) /
        sum(RELATIVE_CONSUMPTION*DRINKERS))

  PC %<>%
    add_column(GAMMA_CONSTANT = sapply(gc[PC$GENDER], `[[`, 1))

  PC %<>%
    mutate(
      GAMMA_SHAPE = 1/GAMMA_CONSTANT,
      GAMMA_SCALE = GAMMA_CONSTANT*PCC_AMONG_DRINKERS)

  PC %<>%
    mutate(
      BB = bb[GENDER][[1]],
      LB = lb,
      UB = ub)

  PC[, "NC"] <- 0
  for(n in 1:nrow(PC)) {
    GAMMA_SHAPE <- PC[[n, "GAMMA_SHAPE"]]
    GAMMA_SCALE <- PC[[n, "GAMMA_SCALE"]]
    LB <- PC[[n, "LB"]]
    UB <- PC[[n, "UB"]]
    imgamma <- function(x) {
      dgamma(x, shape = GAMMA_SHAPE, scale = GAMMA_SCALE)
    }
    PC[[n, "NC"]] <- integrate(imgamma, lower = LB, upper = UB)$value
  }

  PC %<>%
    mutate(DF = P_CD / NC)

  PC[, "P_BAT"] <- 0
  for(n in 1:nrow(PC)) {
    GAMMA_SHAPE <- PC[[n, "GAMMA_SHAPE"]]
    GAMMA_SCALE <- PC[[n, "GAMMA_SCALE"]]
    BB <- PC[[n, "BB"]]
    UB <- PC[[n, "UB"]]
    DF <- PC[[n, "DF"]]
    dfgamma <- function(x) {
      DF * dgamma(x, shape = GAMMA_SHAPE, scale = GAMMA_SCALE)
    }
    PC[[n, "P_BAT"]] <- integrate(dfgamma, lower = BB, upper = UB)$value
  }
  PC %<>%
    mutate(
      R1 = (P_CD - P_BD)  / (P_CD - P_BAT),
      R2 = (P_BD - P_BAT) / (P_CD - P_BAT),
      AAF_FD = 0) %>%
    add_column(
      N_GAMMA = zero,
      INTGRND = zero,
      AAF_CMP = zero)

  for(n in 1:nrow(PC)) {
    PC[[n, "N_GAMMA"]] <- normalized_gamma_factory(PC[n, ])
  }

  PC
}

#' Factory for normalized gamma distributions
#'
#'@param pc_specs is a tibble that contains the variables:
#'  GAMMA_SHAPE: dbl
#'  GAMMA_SCALE: dbl
#'  DF: dbl
#'
#'@return a function object that repreents the normalized gamma distribution
#'  (defined as defined as DF*gamma(x,shape,scale))
#'

normalized_gamma_factory <- function(pc_specs) {
  force(pc_specs)
  GAMMA_SHAPE <- pc_specs[["GAMMA_SHAPE"]]
  GAMMA_SCALE <- pc_specs[["GAMMA_SCALE"]]
  DF <- pc_specs[["DF"]]
  function(x) {
    DF * dgamma(x, shape = GAMMA_SHAPE, scale = GAMMA_SCALE)
  }
}


#' Derives count data for use in calibration and portioning
#'
#'@description
#' -Melts Deaths/Hosps variables into outcome/count (provides better abstraction
#' for calibration methods)
#' -Derives as count/drinkers the region/year/cohort prevalence among drinkers
#' -
#'
#'@param dh Deaths/Hosps as returned by format_v0_dh
#'@param pc Prev/Cons as returned by derive_v0_pc
#'
#'@importFrom magrittr %>% %<>%
#'
#'@export

derive_v0_dh <- function(dh, pc) {
  PC <- pc[
    c("REGION", "YEAR", "GENDER", "AGE_GROUP", "DRINKERS",
      "BB", "LB", "UB", "N_GAMMA")]
  DH <- dh %>%
    dplyr::inner_join(
      PC,
      by = c("REGION", "YEAR", "GENDER", "AGE_GROUP")
    )
  DH
}

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

#### Handle missing data -------------------------------------------------------

#' What value is imputed when a variable is missing?
#'
#' In order to normalize code in future methods, missing data is removed either
#' throug imputation when possible, or by throwing an error when the data is
#' necessary.  impute_with is called when a missing value \emph{can} be imputed.
#'
#' @param var is a string -- a variable name.
#'
#' @return A string or integer that is appropriate for imputation.
#'

impute_with <- function(var) {
  if(grepl("B[0-9]", var)) {var = "BETA"}
  switch(var,
         IM = "(0).(0)",
         CONDITION = "Unspecified",
         REGION = "Unspecified",
         YEAR = "Unspecified",
         AGE_GROUP = "Unspecified",
         OUTCOME = "Combined",
         FUNCTION = "FP",
         RR_FD = 1,
         BINGEF = 1,
         BETA = 0)
}

#' Definition of missing data for the purposes of intermahpr
#'

is.missing <- function(obs) {
  is.na(obs) | is.null(obs) | obs == "."
}

#### Verify and wrangle data ---------------------------------------------------

#' Format relative risk input data to our desired specifications
#'

format_v0_rr <- function(rr) {
  RR <- tibble::as.tibble(rr)
  names(RR) <- do.call(stringr::str_to_upper, list(names(RR)))
  BETAS <- do.call(paste0, list(rep("B", 16), 1:16))
  EXPECTED <- c("IM", "CONDITION", "GENDER", "OUTCOME", "RR_FD", "BINGEF",
                    "FUNCTION", BETAS)
  MISSING <- EXPECTED[!(EXPECTED %in% names(RR))]

  RR[, MISSING] <- NA

  for(var in names(RR)) {
    RR[var][is.missing(RR[var])] <- impute_with(var)
  }

  RR <- readr::type_convert(RR)
}

#' Format prevalence and consumption input data to our desired specifications
#'
#' expected variables:
#'  Needed:
#'    Proportion, Former Drinkers
#'    Proportion, Binge Drinkers
#'    Proportion, Lifetime Abstainers
#'    Proportion, Current Drinkers
#'    Gamma Distribution Parameters:
#'      Per Capita Consumption
#'      Population
#'      Relative Consumption
#'      Correction Factor
#'    Gender (in order to pair with RR curves)
#'  Imputable as "Unspecified":
#'    Region
#'    Year
#'    Gender
#'    Age Group
#'

format_v0_pc <- function(pc) {
  PC <- tibble::as.tibble(pc)
  names(PC) <- do.call(stringr::str_to_upper, list(names(PC)))


  EXPECTED <- c("REGION", "YEAR", "GENDER", "AGE_GROUP", "POPULATION",
                    "PCC_LITRES_YEAR", "CORRECTION_FACTOR",
                    "RELATIVE_CONSUMPTION", "P_LA", "P_FD", "P_CD", "P_BD")
  MISSING <- EXPECTED[!(EXPECTED %in% names(PC))]
  NEEDED <- EXPECTED[c(3,5:12)]
  MISSING_NEEDED <- intersect(MISSING, NEEDED)
  if(length(MISSING_NEEDED) > 0) {
    stop(paste0("Missing and needed variables from the supplied prevalence/",
                "consumption sheet: "),
                paste(MISSING_NEEDED, collapse = ", "))
  }

  PC[, MISSING] <- NA

  for(var in names(PC)) {
    PC[var][is.missing(PC[var])] <- impute_with(var)
  }

  readr::type_convert(PC)
}

#### Process formatted data ----------------------------------------------------

#' Derives relative risk curves from relative risk data
#'
#'@param rr tibble as formatted above
#'@param ext logical indicator of extrapolation, (TRUE = linear, FALSE = capped)
#'
#'@return tibble which is RR with an additional column that contains all
#'  relevent relative risk curves under the variable CURVES. The CURVES variable
#'  has the following structure:
#'    RR[["CURVES"]][[i]] is a list
#'    RR[["CURVES"]][[i]][["BASE_RR"]] is the base relative risk function (i.e.
#'      no extrapolation after 100/150)
#'    RR[["CURVES"]][[i]][["LNXT_RR"]] is the extrapolated relative risk curve
#'    RR[["CURVES"]][[i]][["BNGD_RR"]] is the extrapolated relative risk curve
#'      for binge drinkers

derive_v0_rr <- function(rr, ext) {
  rr[, "EXT"] <- ext
  RR_FN_LIST <- list(apply(rr, 1, function(obs) compile_rr(obs)))
  rr["CURVES"] <- RR_FN_LIST
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
#'@param PC is assumed type tibble with names(PC) = c(YEAR, REGION, GENDER
#' , AGE_GROUP, POPULATION, PCC_LITRES_YEAR, CORRECTION_FACTOR,
#' RELATIVE_CONSUMPTION, P_LA, P_FD, P_CD, P_BD).
#' This function is internally called only, not exported.
#'
#' Note that the bundled data set pc_default satisfies these constraints.
#'
#'@param bb User supplied gender stratified binge barrier
#'@param lb
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
#'@importFrom magrittr "%>%"
#'@importFrom dplyr group_by mutate
#'

derive_v0_pc <- function(pc,
                         bb = c(53.8, 67.25),
                         lb = c(0.03, 0.03),
                         ub = c(250,  250)) {
  PC <- pc %>%
    group_by(REGION, YEAR) %>%
    mutate(PCC_G_DAY = PCC_LITRES_YEAR * 1000 *
             0.7893 * CORRECTION_FACTOR * 0.002739726,
           DRINKERS = POPULATION * P_CD,
           PCAD = PCC_G_DAY * sum(POPULATION) / sum(DRINKERS),
           PCC_AMONG_DRINKERS = RELATIVE_CONSUMPTION * PCAD * sum(DRINKERS) /
             sum(RELATIVE_CONSUMPTION*DRINKERS))
  PC <- PC %>%
    mutate(GAMMA_CONSTANT = ifelse(GENDER == "Female", 1.258^2, 1.171^2),
           GAMMA_SHAPE = 1/GAMMA_CONSTANT,
           GAMMA_SCALE = GAMMA_CONSTANT*PCC_AMONG_DRINKERS,
           INDEX = ifelse(GENDER == "Female", as.integer(1), as.integer(2)))
  PC <- PC %>%
    mutate(BB = bb[INDEX],
           LB = lb[INDEX],
           UB = ub[INDEX])

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

  PC <- PC %>%
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
  PC <- PC %>%
    mutate(R1 = (P_CD - P_BD)  / (P_CD - P_BAT),
           R2 = (P_BD - P_BAT) / (P_CD - P_BAT))
  PC
}

#### Add AAF computing functions -----------------------------------------------

#' Prepares PC and RR data for use in AAF methods
#'
#'@description Performs a full_join of RR and PC over GENDER and adds N_GAMMA
#'  (Normalized GAMMA distribution) function to the CURVES list (N_GAMMA defined
#'  as DF*gamma(x,shape,scale))
#'
#'@param RR relative risk tibble as produced by derive_v*_rr
#'@param PC prevalence and consumption tibble as produced by derive_v*_pc
#'
#'@return tibble with one row per unique region.year.gender.age_group.im combn
#'

join_pc_rr <- function(pc, rr) {
  JOINT <- inner_join(pc, rr, by = "GENDER")
  JOINT[["CURVES"]][[, "N_GAMMA"]] <- function(x) 0
  for(n in 1:nrow(JOINT)) {
    DF <- JOINT[["DF"]][[n]]
    GAMMA_SHAPE <- JOINT[["GAMMA_SHAPE"]][[n]]
    GAMMA_SCALE <- JOINT[["GAMMA_SCALE"]][[n]]
    N_GAMMA <- function(x) {
      DF * dgamma(x, shape = GAMMA_SHAPE, scale = GAMMA_SCALE)
    }
    JOINT[["CURVES"]][[n]][["N_GAMMA"]] <- N_GAMMA
  }
}

#' Obtains AAF computing functions for each observation
#'
#'@param JOINT tibble as produced by join_pc_rr
#'
#'@return tibble with additional variable AAF_FD and CURVES variable wherein
#'  each list contains the additional curve INTGRND and function AAF_CMP.
#'  INTGRND is described in int_factory and AAF_CMP is described in aaf_factory.
#'

add_aaf_fns <- function(joint) {
  joint[["CURVES"]][[, "INTGRND"]] <- function(x) 0
  INTGRND_LIST <- apply(joint, 1, function(obs) int_factory(obs))
  for(n in 1:nrow(joint)) {

  }


  joint[["CURVES"]][[, "AAF_CMP"]] <- function(a,b) 0

  joint[["CURVES"]][[, "INTGRND"]] <- function(x) 0



  for(n in 1:nrow(joint)) {

  }


  joint[["CURVES"]][[, "AAF_CMP"]] <- function(a,b) 0
}


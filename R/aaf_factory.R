#' Factory for normalized gamma distributions
#'
#'@param aaf_specs is a tibble that contains the variables:
#'  GAMMA_SHAPE: dbl
#'  GAMMA_SCALE: dbl
#'  DF: dbl
#'
#'@return a function object that repreents the normalized gamma distribution
#'  (defined as defined as DF*gamma(x,shape,scale))
#'

normalized_gamma_factory <- function(aaf_specs) {
  force(aaf_specs)
  GAMMA_SHAPE <- aaf_specs[["GAMMA_SHAPE"]]
  GAMMA_SCALE <- aaf_specs[["GAMMA_SCALE"]]
  DF <- aaf_specs[["DF"]]
  function(x) {
    DF * dgamma(x, shape = GAMMA_SHAPE, scale = GAMMA_SCALE)
  }
}


#' Factory for aaf integrands
#'
#'@description Given function data (gamma distribution specs, relative risk
#'  curves, binge ratios) integrand_factory produces a function for use as an
#'  integrand in AAF computations
#'
#'@param aaf_specs is a tibble that contains the variables:
#'  BB: dbl
#'  LB: dbl
#'  UB: dbl
#'  R1: dbl
#'  R2: dbl
#'  CURVES:: list of functions:
#'    BASE_RR
#'    LNXT_RR
#'    BNGD_RR
#'    N_GAMMA
#'
#'@return a function object that represents the generalized integrand found in
#'  the generalized AAF computation outlined in the the InterMAHP user guide
#'

integrand_factory <- function(aaf_specs) {
  force(aaf_specs)
  BB <- aaf_specs[["BB"]]
  LB <- aaf_specs[["LB"]]
  UB <- aaf_specs[["UB"]]
  R1 <- aaf_specs[["R1"]]
  R2 <- aaf_specs[["R2"]]
  CURVES <- aaf_specs[["CURVES"]][[1]]
  BASE_RR <- CURVES[["BASE_RR"]]
  LNXT_RR <- CURVES[["LNXT_RR"]]
  BNGD_RR <- CURVES[["BNGD_RR"]]
  N_GAMMA <- CURVES[["N_GAMMA"]]
  NB_INTEGRAND <- function(x) N_GAMMA(x) * (LNXT_RR(x) - 1)
  BD_INTEGRAND <- function(x) N_GAMMA(x) * (BNGD_RR(x) - 1)
  function(x) {
    (x <= BB)*(R1*NB_INTEGRAND(x) + R2*BD_INTEGRAND(x)) +
    (x >  BB)*BD_INTEGRAND(x)
  }
}

#' Factory for alcohol attributable fraction computing functions
#'
#'@description Given function data (gamma distribution specs, relative risk
#'  curves, binge ratios) aaf_factory produces a function f that computes aaf
#'  aafs from a given interval
#'
#'@param aaf_specs is a tibble that contains the variables:
#'  LB: dbl
#'  UB: dbl
#'  RR_FD: dbl
#'  P_FD: dbl
#'  CURVES: list
#'    INTGRND: function
#'
#'@return a function object that takes as input two values a, b such that
#'  LB < a < b < UB and returns as output the alcohol attributable fraction for
#'  outcome OUTCOME due to condition CONDITION among cohort GENDER * AGE_GROUP
#'  in the place and time REGION * YEAR (unused metadata).
#'

aaf_factory <- function(aaf_specs) {
  force(aaf_specs)
  LB <- aaf_specs[["LB"]]
  UB <- aaf_specs[["UB"]]
  RR_FD <- aaf_specs[["RR_FD"]]
  P_FD <- aaf_specs[["P_FD"]]
  INTGRND <- aaf_specs[["CURVES"]][[1]][["INTGRND"]]
  CD_FACTOR <- integrate(INTGRND, lower = LB, upper = UB)$value
  FD_FACTOR <- P_FD * (RR_FD - 1)
  DENOMINATOR <- 1 + FD_FACTOR + CD_FACTOR
  function(a,b) {
    NUMERATOR <- integrate(INTGRND, lower = a, upper = b)$value
    FRACTION <- NUMERATOR / DENOMINATOR
    FRACTION
  }
}

#' Computes alcohol attributable fraction for former drinkers
#'
#'@description Given risk, population, and current drinker integrand information
#'  aaf_fd computes the alcohol attributable fraction among former drinkers
#'
#'@param aaf_specs is a tibble that contains the variables:
#'  RR_FD: dbl
#'  P_FD: dbl
#'  CURVES: list
#'    INTGRND: function
#'
#'@return double to be interpreted as a fraction
#'

aaf_fd <- function(aaf_specs) {
  RR_FD <- aaf_specs[["RR_FD"]]
  P_FD <- aaf_specs[["P_FD"]]
  INTGRND <- aaf_specs[["CURVES"]][[1]][["INTGRND"]]
  CD_FACTOR <- integrate(INTGRND, lower = LB, upper = UB)$value
  FD_FACTOR <- P_FD * (RR_FD - 1)
  DENOMINATOR <- 1 + FD_FACTOR + CD_FACTOR
  FRACTION <- FD_FACTOR / DENOMINATOR
  FRACTION
}

#' Compile a list of AAF functions and values
#'
#'@param aaf_specs is a tibble that contains the variables:
#'  GAMMA_SHAPE: dbl
#'  GAMMA_SCALE: dbl
#'  DF: dbl
#'  BB: dbl
#'  LB: dbl
#'  UB: dbl
#'  R1: dbl
#'  R2: dbl
#'  RR_FD: dbl
#'  P_FD: dbl
#'  CURVES:: list of functions:
#'    BASE_RR
#'    LNXT_RR
#'    BNGD_RR
#'
#'@return a list of the same structure as aaf_specs with no fields altered save
#'  CURVES, which includes the functions N_GAMMA, INTGRND, and AAF_CMP and with
#'  a new variable AAF_FD
#'

compile_aaf <- function(aaf_specs) {
  aaf_specs[["CURVES"]][[1]][["N_GAMMA"]] <- normalized_gamma_factory(aaf_specs)
  print(aaf_specs[["CURVES"]][[1]][["N_GAMMA"]])
  aaf_specs[["CURVES"]][[1]][["INTGRND"]] <- integrand_factory(aaf_specs)
  aaf_specs[["CURVES"]][[1]][["AAF_CMP"]] <- aaf_factory(aaf_specs)
  aaf_specs[["AAF_FD"]] <- aaf_fd(aaf_specs)
}


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
#'  LNXT_RR: fn
#'  BNGD_RR: fn
#'  N_GAMMA: fn
#'
#'@return a function object that represents the generalized integrand found in
#'  the generalized AAF computation outlined in the the InterMAHP user guide
#'

intgrnd_factory <- function(aaf_specs) {
  force(aaf_specs)
  BB <- aaf_specs[["BB"]]
  LB <- aaf_specs[["LB"]]
  UB <- aaf_specs[["UB"]]
  R1 <- aaf_specs[["R1"]]
  R2 <- aaf_specs[["R2"]]
  LNXT_RR <- aaf_specs[["LNXT_RR"]][[1]]
  BNGD_RR <- aaf_specs[["BNGD_RR"]][[1]]
  N_GAMMA <- aaf_specs[["N_GAMMA"]][[1]]
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
#'  INTGRND: fn
#'
#'@return a function object that takes as input a vector x of values >= 0.03 and
#'  returns the alcohol attributable fraction between 0.03 and x for outcome
#'  OUTCOME due to condition CONDITION among cohort GENDER * AGE_GROUP in the
#'  place and time REGION * YEAR
#'

aaf_cmp_factory <- function(aaf_specs) {
  force(aaf_specs)
  LB <- aaf_specs[["LB"]]
  UB <- aaf_specs[["UB"]]
  RR_FD <- aaf_specs[["RR_FD"]]
  P_FD <- aaf_specs[["P_FD"]]
  INTGRND <- aaf_specs[["INTGRND"]][[1]]
  CD_FACTOR <- integrate(INTGRND, lower = LB, upper = UB)$value
  FD_FACTOR <- P_FD * (RR_FD - 1)
  DENOMINATOR <- 1 + FD_FACTOR + CD_FACTOR
  function(x) {
    integral_up_to <- function(up_to) {
      integrate(INTGRND, lower = LB, upper = up_to)$value
    }
    NUMERATOR <- sapply(x, integral_up_to)
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
#'  LB: dbl
#'  UB: dbl
#'  P_FD: dbl
#'  RR_FD: dbl
#'  INTGRND: fn
#'
#'@return double to be interpreted as a fraction
#'

aaf_fd <- function(aaf_specs) {
  LB <- aaf_specs[["LB"]]
  UB <- aaf_specs[["UB"]]
  P_FD <- aaf_specs[["P_FD"]]
  RR_FD <- aaf_specs[["RR_FD"]]
  INTGRND <- aaf_specs[["INTGRND"]][[1]]
  CD_FACTOR <- integrate(INTGRND, lower = LB, upper = UB)$value
  FD_FACTOR <- P_FD * (RR_FD - 1)
  DENOMINATOR <- 1 + FD_FACTOR + CD_FACTOR
  FRACTION <- FD_FACTOR / DENOMINATOR
  FRACTION
}


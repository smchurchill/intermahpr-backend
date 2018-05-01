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
    NUMERATOR <- vapply(x, integral_up_to, 0)
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

#### 100% Calibration ----------------------------------------------------------


#'
#'
#'@param aaf_specs is a tibble that contains the variables:
#'  AAF_CMP: fn
#'  AAF_TOTAL: dbl
#'  P_FD: dbl
#'  RR_FD: dbl
#'  INTGRND: fn
#'
#'

scaled_aaf_cmp_factory <- function(aaf_specs) {
  AAF_CMP <- aaf_specs[["AAF_CMP"]][[1]]
  RESCALE <- 1/aaf_specs[["AAF_TOTAL"]]
  function(x) {
    RESCALE * AAF_CMP(x)
  }
}

#'
#'@param aaf_specs is a tibble that contains the variables:
#'  AAF_CMP: fn
#'

copy_aaf_cmp_factory <- function(aaf_specs) {
  return(aaf_specs[["AAF_CMP"]][[1]])
}


#' Produce a Conditional Probability curve from the given input
#'
#'@description
#'Given a list that contains the true count among a given cohort of
#'a wholly alocohol-attributable condition, the number of drinkers in that
#'cohort, the conditions's IM code, and the gamma function detailing consumption
#'among the cohort, this function calibrates a slope for a loglinear function
#'that estimates the conditional probablity of developing the given condition.
#'The function should be interpreted as conditional probability mass.  I.e., the
#'integral of the function over the relevant range is the probability that a
#'drinker will be afflicted by an event of the given condition over the period
#'of time that the given Count variable was collected.
#'
#'Uses a nonlinear optimizer (COBYLA) to find a loglinear slope for the
#'function f(x) = max(1, exp(k(x-t))) that mimizes the difference between
#'integral(N_GAMMA * (f-1), 0.03, UB) and yearly prevalence (count/drinkers).
#'
#'The goal is to produce a continuous analogue to the relative risk curve for
#'conditions that are wholly attributable to alcohol. The assumption is made
#'that such a condition has a loglinear thresholded (i.e. f(x)=1 for x<t)
#'conditional probablity function on the interval of concern (0.03 to UB grams
#'of ethanol/day, averaged over 1yr).
#'
#'This conditional probability is used to portion a 1.00 AAF_TOTAL among the
#'drinking population.
#'
#'@param rr_specs is a tibble that contains the variables:
#'  IM: chr
#'  COUNT: dbl
#'  DRINKERS: dbl
#'  N_GAMMA: fn
#'  LB: dbl
#'  BB: dbl
#'  UB: dbl
#'
#'@return a function whose values are binge modified.  This affects conditions
#'  5.2 and 5.5, ischaemic heart disease and stroke respectively, by removing a
#'  protective effect at low levels of consumption
#'

calibration_factory <- function(dh_specs) {
  ## Determine whether we want to return early.
  IM <- dh_specs[["IM"]]
  COUNT <- dh_specs[["COUNT"]]
  LB <- dh_specs[["LB"]]
  BB <- dh_specs[["BB"]]
  UB <- dh_specs[["UB"]]
  THRESHOLD <- BB
  if(grepl("6", IM)) {
    THRESHOLD <- LB
  }

  ## If there's no count to calibrate against, we return early with a function
  ## that will at least evaluate 1 at the upper bound and 0 below the threshold.
  if(COUNT <= 0) {
    return(
      function(x) {
        (x > THRESHOLD)*(x-THRESHOLD)/(UB - THRESHOLD)
      }
    )
  }

  DRINKERS <- dh_specs[["DRINKERS"]]
  N_GAMMA <- dh_specs[["N_GAMMA"]][[1]]

  MASS <- function(k) {
    function(x) {
      N_GAMMA(x) * (exp(pmax(0, k*(x-THRESHOLD)))-1)
    }
  }

  EST_COUNT <- function(k) {
    DRINKERS*integrate(MASS(k), LB, UB)$value
  }

  EST_ERROR <- function(k) {
    abs(EST_COUNT(k) - COUNT)
  }

  OPTR <- nloptr::nloptr(
    x0 = 0.01,
    eval_f = EST_ERROR,
    lb = 0,
    ub = 1,
    opts = list(
      "algorithm" = "NLOPT_LN_COBYLA",
      "xtol_rel"  = 1.0e-10
    )
  )

  SLOPE <- OPTR$solution

  function(x) {
    INTGRND <- MASS(SLOPE)
    integral_up_to <- function(up_to) {
      DRINKERS*integrate(MASS(SLOPE), lower = LB, upper = up_to)$value/COUNT
    }
    vapply(x, integral_up_to, 0)
  }
}

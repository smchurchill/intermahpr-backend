##### g02-factories ############################################################
##
## Factory functions called from functions in s02-process
##
##
##
##


#### Calibrate AAF Computer Factory --------------------------------------------

#' Slope Factory
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
#'The nonlinear optimization always estimates conservatively.  For conditions
#'with observed incidence at leats 1.8 per 10,000 drinkers, the AAF will round
#'to 1 within 2 decimal places.  There is no such guarantee for conditions with
#'incidence less than 1.8 per 10,000 drinkers.
#'
#'@param IM InterMAHP condition code
#'@param COUNT Count of condition events over the relevant time period
#'@param DRINKERS Estimated # of drinkers in population
#'@param N_GAMMA Normalized gamma distribution (used as exposure mass function)
#'@param LB,BB,UB lower, binge, and upper bounds of consumption
#'
#'@return slope of loglinear conditional probability mass function for risk as a
#'result of exposure
#'
#'@export

slope_calibration <- function(IM, COUNT, DRINKERS, N_GAMMA, LB, BB, UB) {
  ## If there's no count to calibrate against, we return early with a function
  ## that will at least evaluate 1 at the upper bound and 0 below the threshold.
  if(COUNT <= 0) return(0)

  ## Threshold is binge for all conditions except digestive
  THRESHOLD <- BB
  if(grepl("6", IM)) {
    THRESHOLD <- LB
  }

  MASS <- function(k) {
    function(x) {
      N_GAMMA(x) * (exp(pmax(0, k*(x-THRESHOLD)))-1)
    }
  }

  EST_COUNT <- function(k) {
    integrate(MASS(k), LB, UB)$value
  }

  EST_ERROR <- function(k) {
    abs(EST_COUNT(k) - (COUNT/DRINKERS))
  }

  nloptr::nloptr(
    x0 = 0.01,
    eval_f = EST_ERROR,
    lb = 0,
    ub = 1,
    opts = list(
      "algorithm" = "NLOPT_LN_COBYLA",
      "xtol_rel"  = 1.0e-20
    )
  )$solution
}


#'Conditional Probability Factory
#'@description Invokes slope_calibration and wraps the result in a loglinear
#'probability mass function
#'@inheritParams slope_calibration
#'@return Conditional probability mass function for risk incurred as a result of
#'exposure
#'@export

cond_prob_factory <- function(IM, COUNT, DRINKERS, N_GAMMA, LB, BB, UB) {
  slope <- slope_calibration(IM, COUNT, DRINKERS, N_GAMMA, LB, BB, UB)
  ## Threshold is binge for all conditions except digestive
  THRESHOLD <- BB
  if(grepl("6", IM)) {
    THRESHOLD <- LB
  }
  function(x) exp(pmax(0, slope*(x-THRESHOLD)))-1
}

#'Calibrated AAF Factory
#'@description Invokes cond_prob_factory and wraps the result multiplied by the
#'given exposure mass function N_GAMMA
#'@inheritParams slope_calibration
#'@return AAF_CMP function for distributing harm for wholly attributable
#'functions
#'@export

aaf_calibration_factory <- function(IM, COUNT, DRINKERS, N_GAMMA, LB, BB, UB) {
  ## Threshold is binge for all conditions except digestive
  THRESHOLD <- BB
  if(grepl("6", IM)) {
    THRESHOLD <- LB
  }

  ## If there's no count to calibrate against, we return early with a function
  ## that will at least evaluate 1 at the upper bound and 0 below the threshold.
  if(COUNT <= 0) {
    return(function(x) (x > THRESHOLD)*(x-THRESHOLD)/(UB - THRESHOLD))
  }

  COND_PROB <- cond_prob_factory(IM, COUNT, DRINKERS, N_GAMMA, LB, BB, UB)
  MASS <- function(x) DRINKERS*(N_GAMMA %prod% COND_PROB)(x)/COUNT
  integral_up_to <- function(up_to) {
    if(up_to >= UB) return(1)
    if(up_to <= THRESHOLD) return(0)
    integrate(MASS, lower = LB, upper = up_to)$value
  }

  function(x) vapply(x, integral_up_to, 0)
}

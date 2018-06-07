#' Slope Calibration
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
#'@param target Observed incidence to calibrate against
#'@param mass population exposure mass function to calibrate against
#'@param lb lower bound of consumption at which condition occurs
#'@param ub upper bound of consumption

#'@return slope of loglinear conditional probability mass function for risk as a
#'result of exposure
#'
#'@export

calibrateSlope <- function(target, mass, lb, ub) {
  if(target <= 0) return(0)

  integrand <- function(k) function(x) mass(x) * (exp(pmax(0, k*(x-lb)))-1)
  estimate <- function(k) integrate(integrand, lb, ub)$value
  error <- function(k) abs(estimate(k) - target)

  nloptr::nloptr(
    x0 = 0.01,
    eval_f = error,
    lb = 0,
    ub = 1,
    opts = list(
      "algorithm" = "NLOPTR_LN_COBYLA",
      "xtol_rel" = 1.0e-20
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

makeConditionalProbability <- function(slope, lb) {
  function(x) exp(pmax(0, slope*(x-lb)))-1
}

makeCurrentCalibratedFactory <- function(target, mass, lb, ub) {
  slope <- slope_calibration(target = target, mass = mass, lb = lb, ub = ub)
  cond_prob <- cond_prob_factory(slope = slope, lb = lb)

  ## *d_aaf_fct functions have positional ... arguments. positions are:
  ## ..1: population exposure mass function
  ## ..2: proportion of drinkers below BB that do not binge
  ## ..3: proportion of drinkers below BB that do binge
  ## ..4: lower bound
  ## ..5: binge barrier
  ## ..6: upper bound
  ## for dependent aafs, only ..1 is used.
  function(...) {
    cdf_fct_fct(f = ..1 %prod% cprob, lb = lb)
  }
}

makeFormerCalibratedFactory <- function(...) function(...) function(x) 0


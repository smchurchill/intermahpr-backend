#### Incidence Calibraition ----------------------------------------------------

#' Calibrate a loglinear Slope
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
#'
#'@return slope of loglinear conditional probability mass function for risk as a
#'result of exposure
#'

calibrateSlope <- function(target, mass, lb, ub) {
  if(target <= 0) return(0)

  integrand <- function(k) function(x) mass(x) * (exp(pmax(0, k*(x-lb)))-1)
  estimate <- function(k) integrate(integrand(k), lb, ub)$value
  error <- function(k) abs(estimate(k) - target)

  nloptr::nloptr(
    x0 = 0.01,
    eval_f = error,
    lb = 0,
    ub = 1,
    opts = list(
      "algorithm" = "NLOPT_LN_COBYLA",
      "xtol_rel" = 1.0e-20
    )
  )$solution
}


#' Factory for loglinear conditional probability functions
#'@description Invokes slope_calibration and wraps the result in a loglinear
#'probability mass function
#'@return Conditional probability mass function for risk incurred as a result of
#'exposure

makeConditionalProbability <- function(slope, lb) {
  function(x) {
    exp(pmax(0, slope*(x-lb)))-1
  }
}

#' Factory for current drinker's AAF computer factory for a condition with
#'  calibrated incidence estimator
#'

makeCurrentCalibratedFactory <- function(target, clbr_mass, lb, ub) {
  slope <- calibrateSlope(target = target, mass = clbr_mass, lb = lb, ub = ub)
  c_prob <- makeConditionalProbability(slope = slope, lb = lb)
  if(target > 0) reciprocal_target <- 1/target
  else reciprocal_target <- 1
  function(args) {
    integrand <- function(x) reciprocal_target * (args$mass %prod% c_prob)(x)
    makeIntegrator(f = integrand, lb = lb, ub = ub)
  }
}

#' Factory for former drinker's AAF computer factory for a condition with
#'  calibrated incidence estimator
#'
#' The relative risk for former drinkers is 1 by definition, so these functions
#'  are identically zero.
#'

makeFormerCalibratedFactory <- function(...) {
  function(...) {
    function(x) {
      0
    }
  }
}


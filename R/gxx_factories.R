#### Relative Risk Choice & Mod Factories --------------------------------------

#' Choose which relative risk function is appropriate from the given input
#'
#'@description Given a list that contains the string variables IM, GENDER, and
#'  FUNCTION, and sixteen double variables B1 through B16, returns a relative
#'  risk function.
#'  Currently supported levels of IM can be found in the User Guide. Supported
#'  levels for GENDER are Female and Male, and for FUNCTION are FP, Spline, and
#'  Step.
#'  IM, FUNCTION, and GENDER are used to determine which function to produce in
#'  extraordinary cases.  The most common case, however, is the FP (Fractional
#'  Polynomial) model. When FUNCTION = FP and the given condition is not
#'  extraordinary, a relative risk function will still be produced even if
#'  Gender/IM are imputed by intermahpr.
#'
#'@param IM string, Intermahp condition code
#'@param GENDER string
#'@param FUNCTION string, "FP", "Step", or "Spline"
#'@param BETAS dbl, vector of length 16
#'
#'@return a function object that describes a base relative risk curve
#'

base_rr_factory <- function(IM, GENDER, FUNCTION, BETAS) {
  FN_RR <- function(x) {0}

  if(FUNCTION == "FP"){
    FN_RR <- fractional_polynomial_rr(BETAS)
  }

  if(FUNCTION == "Step" & GENDER == "Female") {
    FN_RR <- hiv_f_rr
  }
  if(FUNCTION == "Step" & GENDER == "Male") {
    FN_RR <- hiv_m_rr
  }

  if(FUNCTION == "Spline") {
    if(IM == "(6).(3)") {
      FN_RR <- acute_pancreatitis_f_rr
    }
    else if(GENDER == "Female") {
      FN_RR <- hypertension_f_rr
    }
    else {
      FN_RR <- hypertension_m_rr
    }
  }
  FN_RR
}

#' Produce an extrapolation-modified relative risk curve
#'
#'@description given a relative risk curve and some metadata, produces a
#'  risk curve that has the desired extrapolation behaviour beyond 150 g/day,
#'  or 100 g/day for IHD.
#'
#'@param BASE_RR continuous function in x as returned by base_rr_factory
#'@param X2 extrapolate after x = X2
#'@param Y2 Y2 = BASE_RR(X2)
#'@param SLOPE slope of the extrapolation
#'
#'
#'@return a function whose values above 150 (100 for IHD) are linearly
#'  extrapolated (see InterMAHP guide for details)
#'

linear_extrapolation_factory <- function(BASE_RR, X2, Y2, SLOPE) {
  LINE <- function(x) {
    Y2+ SLOPE*(x-X2)
  }

  function(x) {
    ((0 < x) & (x < X2))*BASE_RR(x) + (X2 <= x)*LINE(x)
  }
}

#' Produce a Relative-Risk-For-Bingers curve from the given input
#'
#'@description given a list that contains the string variable IM and the double
#'  variable BINGEF, produces a Relative-Risk-For-Bingers curve.  IM is used for
#'  the curves for IHD and ischaemic stroke, where any protective J-shape is
#'  forfeited by bingers.  All curves are then multiplied by BINGEF which has
#'  the value of 1 in all cases except for condition classes 7, 8, 9, which
#'  account for motor vehicle collisions, unintentional injuries, and
#'  intentional injuries.
#'
#'@param IM string, Intermahp condition code
#'@param BINGEF double, rescales injuries
#'@param LNXT_RR continuous vector valued function as returned by lin_ext_fact
#'
#'
#'@return a function whose values are binge modified.  This affects conditions
#'  5.2 and 5.5, ischaemic heart disease and stroke respectively, by removing a
#'  protective effect at low levels of consumption
#'

binge_risk_factory <- function(IM, BINGEF, LNXT_RR) {
  if(IM == "(5).(2)" | IM == "(5).(5)") {
    return(function(x) BINGEF*pmax(1, LNXT_RR(x)))
  }
  else {
    return(function(x) BINGEF*LNXT_RR(x))
  }
}

#### Gamma Function Factories --------------------------------------------------

#' Factory for gamma distributions
#'
#'@description
#'We want easier access to gamma distributions with a wide array of fixed params
#'
#'@param GAMMA_SHAPE,GAMMA_SCALE Gamma parameters to be supplied to dgamma
#'
#'

gamma_factory <- function(GAMMA_SHAPE, GAMMA_SCALE) {
  function(x) dgamma(x, shape = GAMMA_SHAPE, scale = GAMMA_SCALE)
}

#' Factory for normalized gamma distributions
#'
#'@param GAMMA Input gamma distribution
#'@param FACTOR Factor to scale GAMMA by
#'
#'

normalized_gamma_factory <- function(GAMMA_SHAPE, GAMMA_SCALE, DF) {
  function(x) DF * dgamma(x, shape = GAMMA_SHAPE, scale = GAMMA_SCALE)
}


#### Base AAF Computer Factories -----------------------------------------------

#' Factory for aaf integrands
#'
#'@description Given function data (gamma distribution specs, relative risk
#'  curves, binge ratios) integrand_factory produces a function for use as an
#'  integrand in AAF computations
#'
#'@param BB dlb, binge barrier
#'@param LB dbl, lower bound of consumption
#'@param UB dbl, upper bound of consumption
#'@param R1 dbl, proportion of drinkers below BB that do not binge
#'@param R2 dbl, proportion of drinkers below BB that do binge
#'@param LNXT_RR fn, extrapolated relative risk for nonbingers
#'@param BNGD_RR fn, extrapolated relative risk for bingers
#'@param N_GAMMA fn, normalized distribution of drinkers
#'
#'@return a function object that represents the generalized integrand found in
#'  the generalized AAF computation outlined in the the InterMAHP user guide
#'

intgrnd_factory <- function(BB, LB, UB, R1, R2, LNXT_RR, BNGD_RR, N_GAMMA) {
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
#'@param LB dbl, lower bound of consumption
#'@param UB dbl, upper bound of consumption
#'@param RR_FD dbl, relative risk for former drinkers
#'@param P_FD dbl, proportion of former drinkers
#'@param INTGRND fn, continuous function to be integrated
#'
#'
#'@return a function object that takes as input a vector x of values >= 0.03 and
#'  returns the alcohol attributable fraction between 0.03 and x for outcome
#'  OUTCOME due to condition CONDITION among cohort GENDER * AGE_GROUP in the
#'  place and time REGION * YEAR
#'

aaf_cmp_factory <- function(LB, UB, RR_FD, P_FD, INTGRND) {
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
#'@param LB dbl, lower bound of consumption
#'@param UB dbl, upper bound of consumption
#'@param RR_FD dbl, relative risk for former drinkers
#'@param P_FD dbl, proportion of former drinkers
#'@param INTGRND fn, continuous function to be integrated
#'
#'@return dbl<1 to be interpreted as a fraction
#'

aaf_fd <- function(LB, UB, RR_FD, P_FD, INTGRND) {
  if(sum(is.na(c(LB, UB, RR_FD, P_FD)))) { return (0.0) }
  CD_FACTOR <- integrate(INTGRND, lower = LB, upper = UB)$value
  FD_FACTOR <- P_FD * (RR_FD - 1)
  DENOMINATOR <- 1 + FD_FACTOR + CD_FACTOR
  FRACTION <- FD_FACTOR / DENOMINATOR
  FRACTION
}

#### Calibrate AAF Computer Factory --------------------------------------------

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

calibration_factory <- function(IM, COUNT, DRINKERS, N_GAMMA, LB, BB, UB) {
  ## Determine whether we want to return early.
  # IM <- dh_specs[["IM"]]
  # COUNT <- dh_specs[["COUNT"]]
  # LB <- dh_specs[["LB"]]
  # BB <- dh_specs[["BB"]]
  # UB <- dh_specs[["UB"]]
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

  # DRINKERS <- dh_specs[["DRINKERS"]]
  # N_GAMMA <- dh_specs[["N_GAMMA"]][[1]]

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


#### Relative Risk Curves ------------------------------------------------------

#' Special Female HIV Relative Risk Function
#'
#'@param x A vector of x values at which to evaluate the RR function
#'


hiv_f_rr <- function(x) {
  ((0 < x) & (x < 49))*1 + (x >= 49)*1.54
}

#' Special Male HIV Relative Risk Function
#'
#'@param x A vector of x values at which to evaluate the RR function
#'


hiv_m_rr <- function(x) {
  ((0 < x) & (x < 61))*1 + (x >= 61)*1.54
}

#' Special Spline Female Hypertension Relative Risk Function
#'
#'@param x A vector of x values at which to evaluate the RR function
#'
#'
#' 0.0025 = 1/400
#'

hypertension_f_rr <- function(x){
  SPLINE =
    (x < 75)*(
      (x > 0)*0*x +
        (x > 18.9517)*(
          -0.0154196*x +
            0.0217586*(
              x^3 + 1*(x-20)^3 - 2*(x-10)^3
            )*0.0025
        )
    )+
    (x >= 75)*(0.9649937)
  exp(SPLINE)
}

#' Special Spline Male Hyptertension Relative Risk Function
#'
#'@param x A vector of x values at which to evaluate the RR function
#'
#' 0.01851852 = 1/54
#' 0.0001777778 = 1/(75^2)
#'

hypertension_m_rr <- function(x){
  SPLINE =
    (x > 0)*(0.0150537*x -
               0.0156155*(
                 x^3 +
                   (x >= 75)*(21*(x-75)^3)*0.01851852 -
                   (x >= 21)*(75*(x-21)^3)*0.01851852
               )*0.0001777778
    )
  exp(SPLINE)
}

#' Special Spline Female Acute Pancreatitis Relative Risk Function
#'
#'@param x A vector of x values at which to evaluate the RR function
#'
#'
#' 0.04 = 1/25
#' 0.0007304602 = 1/(37^2)
#'

acute_pancreatitis_f_rr <- function(x){
  spline =
    (x < 108)*(
      (x > 0)*(-0.0272886)*x +
        0.0611466*(
          (x >= 3)*((x-3)^3)+
            (x >= 40)*(12*((x-40)^3)*0.04)-
            (x >= 15)*(37*((x-15)^3)*0.04)
        )*0.0007304602
    )+
    (x >= 108)*(2.327965)
  exp(spline)
}

#' Get Fractional Polynomial Relative Risk Function
#'
#'@param B The numeric vector of Beta values needed to produce a fractional
#'  polynomial
#'

#### Fractional polynomial Factory ---------------------------------------------

fractional_polynomial_rr <- function(betas) {
  force(betas)
  function(x) {
    M = matrix(0,length(x),16)
    sqrtx = sqrt(x)
    logx = log(x)
    rx = 1 / x
    rsqrtx = 1/sqrtx
    if(betas[[1]]!=0) {M[,1] = betas[[1]] *rx*rx       }
    if(betas[[2]]!=0) {M[,2] = betas[[2]] *rx          }
    if(betas[[3]]!=0) {M[,3] = betas[[3]] *rsqrtx      }
    if(betas[[4]]!=0) {M[,4] = betas[[4]] *logx        }
    if(betas[[5]]!=0) {M[,5] = betas[[5]] *sqrtx       }
    if(betas[[6]]!=0) {M[,6] = betas[[6]] *x           }
    if(betas[[7]]!=0) {M[,7] = betas[[7]] *x*x         }
    if(betas[[8]]!=0) {M[,8] = betas[[8]] *x*x*x       }
    if(betas[[9]]!=0) {M[,9] = betas[[9]] *rx*rx*logx  }
    if(betas[[10]]!=0){M[,10]= betas[[10]]*rx*logx     }
    if(betas[[11]]!=0){M[,11]= betas[[11]]*rsqrtx*logx }
    if(betas[[12]]!=0){M[,12]= betas[[12]]*logx*logx   }
    if(betas[[13]]!=0){M[,13]= betas[[13]]*sqrtx*logx  }
    if(betas[[14]]!=0){M[,14]= betas[[14]]*x*logx      }
    if(betas[[15]]!=0){M[,15]= betas[[15]]*x*x*logx    }
    if(betas[[16]]!=0){M[,16]= betas[[16]]*x*x*x*logx  }
    exp(rowSums(M))
  }
}

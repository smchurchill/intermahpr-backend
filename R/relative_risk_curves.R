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



fp_rr <- function(betas) {
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

#### Curve choosers and modifiers ----------------------------------------------

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
#'@param rr_specs is a tibble that contains the variables:
#'  IM: chr
#'  GENDER: chr
#'  FUNCTION: chr
#'  B1-B16: dbl
#'
#'@return a function object that describes a base relative risk curve
#'

set_rr <- function(rr_specs) {
  force(rr_specs)
  IM <- rr_specs[["IM"]]
  GENDER <- rr_specs[["GENDER"]]
  FUNCTION <- rr_specs[["FUNCTION"]]
  B_INDEX <- do.call(paste0, list(rep("B", 16), 1:16))
  BETAS <- as.numeric(rr_specs[B_INDEX])
  FN_RR <- function(x) {0}

  if(FUNCTION == "FP"){
    FN_RR <- fp_rr(BETAS)
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
#'@param rr_specs is a tibble that contains the variables:
#'  IM: chr
#'  GENDER: chr
#'  EXT: lgl
#'
#'@return a function whose values above 150 (100 for IHD) are linearly
#'  extrapolated (see InterMAHP guide for details)
#'

ext_rr <- function(rr_specs) {
  force(rr_specs)
  IM <- rr_specs[["IM"]]
  GENDER <- rr_specs[["GENDER"]]
  EXT <- rr_specs[["EXT"]]
  FN_RR <- rr_specs[["BASE_RR"]][[1]]

  if(GENDER == "Female" & (IM %in% c("(5).(1)", "(6).(3)"))) {
    EXT <- FALSE
  }

  X1 <- 100 + ifelse(IM == "(5).(2)", -50, 0)
  X2 <- 150 + ifelse(IM == "(5).(2)", -50, 0)
  Y1 <- FN_RR(X1)
  Y2 <- FN_RR(X2)
  SLOPE <- ifelse(EXT, (Y2-Y1)/(X2-X1), 0)

  LINE <- function(x) {
    Y2+ SLOPE*(x-X2)
  }

  function(x) {
    ((0 < x) & (x < X2))*FN_RR(x) + (X2 <= x)*LINE(x)
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
#'@param rr_specs is a tibble that contains the variables:
#'  IM: chr
#'  BINGEF: dbl
#'
#'@return a function whose values are binge modified.  This affects conditions
#'  5.2 and 5.5, ischaemic heart disease and stroke respectively, by removing a
#'  protective effect at low levels of consumption
#'

bng_rr <- function(rr_specs) {
  force(rr_specs)
  IM <- rr_specs[["IM"]]
  BINGEF <- rr_specs[["BINGEF"]]
  TMP_RR <- rr_specs[["LNXT_RR"]][[1]]
  FNC_RR <- rr_specs[["LNXT_RR"]][[1]]
  if(IM %in% c("(5).(2)","(5).(5)")) {
    FNC_RR <- function(x) pmax(1, TMP_RR(x))
  }
  function(x) BINGEF*FNC_RR(x)
}

#### 100% Calibration ----------------------------------------------------------

#' Slope Solver for Thresholded Loglinear Functions
#'
#' Uses a nonlinear optimizer (COBYLA) to find a loglinear slope for the
#' function f(x) = max(1, exp(k(x-t))) that mimizes the difference between
#' integral(N_GAMMA * (f-1), 0.03, UB) and yearly prevalence (count/drinkers).
#'
#' The goal is to produce a continuous analogue to the relative risk curve for
#' conditions that are wholly attributable to alcohol. The assumption is made
#' that such a condition has a loglinear thresholded (i.e. f(x)=1 for x<t)
#' conditional probablity function on the interval of concern (0.03 to UB grams
#' of ethanol/day, averaged over 1yr).
#'
#' This conditional probability is used to portion a 1.00 AAF_TOTAL among the
#' drinking population.
#'
#'@param pc_dh_specs is a tibble that contains the variables:
#'  IM: chr
#'  COUNT: dbl
#'  N_GAMMA: fn
#'  LB: dbl
#'  BB: dbl
#'  UB: dbl
#'  DRINKERS: dbl
#'

fit_slope <- function(pc_dh_specs) {
  IM <- pc_dh_specs[["IM"]]
  COUNT <- pc_dh_specs[["COUNT"]]
  N_GAMMA <- pc_dh_specs[["N_GAMMA"]]
  LB <- pc_dh_specs[["LB"]]
  BB <- pc_dh_specs[["BB"]]
  UB <- pc_dh_specs[["UB"]]
  DRINKERS <- pc_dh_specs[["DRINKERS"]]

  THRESHOLD <- BB
  if(grepl("6", IM)) {
    THRESHOLD <- LB
  }

  FN <- function(k) {
    function(x) {
      N_GAMMA(x) * (exp(pmax(0, k*(x-THRESHOLD))) - 1)
    }
  }

  EST_COUNT <- function(k) {

  }

}

#### Relative Risk Choice & Mod Factories --------------------------------------

#' Factory for base well-defined rel. risk
#'
#' Choose which relative risk function is appropriate from the given input
#'
#'@param im string, Intermahp condition code
#'@param gender string
#'@param form string, "FP", "Step", or "Spline"
#'@param betas dbl, vector of length 16
#'
#'@return a function object that describes a base relative risk curve
#'

makeBaseRisk <- function(im, gender, form, betas) {
  if(form == "FP") return(makeFractionalPolynomial(betas))
  if(form == "Step") return(makeHivRisk(gender))
  if(im == "(6).(3)" & form == "Spline") return(makeAcutePancreatitisRisk())
  if(im == "(5).(1)" & form == "Spline") return(makeHypertensionRisk(gender))
  return(function(x) 1)
}

#' Factory for linearly extrapolated well-defined rel. risk
#'
#'@description
#'  Produces an extrapolation-modified relative risk curve.
#'  Given a relative risk curve and some metadata, produces a
#'  risk curve that has the desired extrapolation behaviour beyond 150 g/day,
#'  or 100 g/day for IHD.
#'
#'@param base_risk continuous function in x as returned by base_rr_factory
#'@param x2 extrapolate after x = X2
#'@param y2 Y2 = BASE_RR(X2)
#'@param slope slope of the extrapolation
#'
#'
#'@return a function whose values after x2 are are linearly extrapolated
#'  (see InterMAHP guide for details)
#'

makeExtrapolatedRisk <- function(base_risk, x2, y2, slope) {
  line <- function(x) y2+ slope*(x-x2)

  function(x) {
    ((0 < x) & (x < x2))*base_risk(x) + (x2 <= x)*line(x)
  }
}

#' Factory for binge-modified well-defined rel. risk.
#'
#'@description
#' Produces a Relative-Risk-For-Bingers curve from the given input
#' Given a list that contains the string variable IM and the double
#'  variable BINGEF, produces a Relative-Risk-For-Bingers curve.  IM is used for
#'  the curves for IHD and ischaemic stroke, where any protective J-shape is
#'  forfeited by bingers.  All curves are then multiplied by BINGEF which has
#'  the value of 1 in all cases except for condition classes 7, 8, 9, which
#'  account for motor vehicle collisions, unintentional injuries, and
#'  intentional injuries.
#'
#'@param im string, Intermahp condition code
#'@param bingef double, rescales injuries
#'@param ext_risk continuous vector valued function as returned by lin_ext_fact
#'
#'
#'@return a function whose values are binge modified.  This affects conditions
#'  5.2 and 5.5, ischaemic heart disease and stroke respectively, by removing a
#'  protective effect at low levels of consumption
#'

makeBingeRisk <- function(im, bingef, ext_risk) {
  min_risk <- 0
  if(im == "(5).(2)" | im == "(5).(5)") min_risk <- 1

  function(x) bingef*pmax(min_risk, ext_risk(x))
}


#### Special Relative Risk Curves ----------------------------------------------

#' Factory for Special HIV Relative Risk Functions
#'
#'@param gender decides which function to return
#'

makeHivRisk  <- function(gender) {
  if(gender == "Female") {
    return(function(x) ((0 < x) & (x < 49))*1 + (x >= 49)*1.54)
  }
  else if(gender == "Male") {
    return(function(x) ((0 < x) & (x < 61))*1 + (x >= 61)*1.54)
  }
  else {
    stop("The following gender has no hardcoded HIV risk: ", gender)
  }
}


#' The spline for female hypertension
#'
#' 0.0025 = 1/400
#'

femaleHypertensionSpline <- function(x) {
  (x < 75)*(
    (x > 0)*0*x +
      (x > 18.9517)*(
        -0.0154196*x +
          0.0217586*(
            x^3 + 1*(x-20)^3 - 2*(x-10)^3
          )*0.0025
      )
  )+(x >= 75)*(0.9649937)
}

#' The spline for male hypertension
#'
#' 0.04 = 1/25
#' 0.0007304602 = 1/(37^2)
#'

maleHypertensionSpline <- function(x) {
  (x > 0)*(
    0.0150537*x -
      0.0156155*(
        x^3 +
          (x >= 75)*(21*(x-75)^3)*0.01851852 -
          (x >= 21)*(75*(x-21)^3)*0.01851852
      )*0.0001777778
  )
}

#' Factory for Special Spline Hypertension Relative Risk Functions
#'
#'@param gender decides which function to return
#'

makeHypertensionRisk <- function(gender) {
  if(gender == "Female") {
    return(function(x) exp(femaleHypertensionSpline(x)))
  }
  else if (gender == "Male") {
    return(function(x) exp(maleHypertensionSpline(x)))
  }
  else{
    stop("The following gender has no hardcoded Hypertension risk: ", gender)
  }

}

#' The spline for Female Acute Pancreatitis
#'
#' 0.04 = 1/25
#' 0.0007304602 = 1/(37^2)
#'

FemaleAcutePancreatitisSpline <- function(x){
  (x < 108)*(
    (x > 0)*(-0.0272886)*x +
      0.0611466*(
        (x >= 3)*((x-3)^3)+
          (x >= 40)*(12*((x-40)^3)*0.04)-
          (x >= 15)*(37*((x-15)^3)*0.04)
      )*0.0007304602
  )+
    (x >= 108)*(2.327965)
}

#' Factory for Special Spline Female Acute Pancreatitis Relative Risk Functions
#'

makeAcutePancreatitisRisk <- function() {
  return(function(x) exp(FemaleAcutePancreatitisSpline(x)))
}

#### Fractional Polynomial Relative Risk Curves --------------------------------

#' The master list of fractional polynomials
#'

fractional_polynomials <- list(
    function(x) 1 / x / x,
    function(x) 1 / x,
    function(x) 1 / sqrt(x),
    function(x) log(x),
    function(x) sqrt(x),
    function(x) x,
    function(x) x*x,
    function(x) x*x*x,
    function(x) log(x) / x / x,
    function(x) log(x) / x,
    function(x) log(x) / sqrt(x),
    function(x) log(x)^2,
    function(x) log(x) * sqrt(x),
    function(x) x*log(x),
    function(x) x*x*log(x),
    function(x) x*x*x*log(x)
)


#' Factory for Fractional Polynomial Relative Risk Functions
#'
#'@param betas The numeric vector of Beta values needed to produce a fractional
#'  polynomial
#'

makeFractionalPolynomial <- function(betas) {

  nonzero_fp <- fractional_polynomials[betas != 0]
  nonzero_b <- betas[betas != 0]

  if(length(nonzero_b) == 0) {return(function(...) 1)}

  function(x) {
    exp(
      Reduce(
        x = lapply(
          X = nonzero_fp,
          FUN = function(f) f(x)
        ),
        f = cbind,
        init = numeric(0)
      ) %*%
        nonzero_b
    )
  }
}

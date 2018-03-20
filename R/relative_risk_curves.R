#' Special Female HIV Relative Risk Function
#'
#'@param x A vector of x values at which to evaluate the RR function
#'


hiv_f_rr <- function(x) {
  ((0<x)&(x<49))*1 + (x>=49)*1.54
}

#' Special Male HIV Relative Risk Function
#'
#'@param x A vector of x values at which to evaluate the RR function
#'


hiv_m_rr <- function(x) {
  ((0<x)&(x<61))*1 + (x>=61)*1.54
}

#' Special Spline Female Hypertension Relative Risk Function
#'
#'@param x A vector of x values at which to evaluate the RR function
#'
#'
#' 0.0025 = 1/400
#'

hypertension_f_rr <- function(x){

  spline =
    (x<75)*(
      (x> 0)*0*x +
        (x>18.9517)*(
          -0.0154196*x +
            0.0217586*(
              x^3 + 1*(x-20)^3 - 2*(x-10)^3
            )*0.0025
        )
    )+
    (x>=75)*(0.9649937)
  exp(spline)
}

#' Special Spline Male Hyptertension Relative Risk Function
#'
#'@param x A vector of x values at which to evaluate the RR function
#'
#' 0.01851852 = 1/54
#' 0.0001777778 = 1/(75^2)
#'

hypertension_m_rr <- function(x){
  spline =
    (x>0)*(0.0150537*x -
             0.0156155*(
               x^3 +
                 (x>=75)*(21*(x-75)^3)*0.01851852 -
                 (x>=21)*(75*(x-21)^3)*0.01851852
             )*0.0001777778
    )
  exp(spline)
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
    (x<108)*(
      (x>0)*(-0.0272886)*x +
        0.0611466*(
          (x>= 3)*((x-3)^3)+
            (x>=40)*(12*((x-40)^3)*0.04)-
            (x>=15)*(37*((x-15)^3)*0.04)
        )*0.0007304602
    )+
    (x>=108)*(2.327965)
  exp(spline)
}

#' Get Fractional Polynomial Relative Risk Function
#'
#'@param B The numeric vector of Beta values needed to produce a fractional
#'  polynomial
#'



fp_rr <- function(B) {
  force(B)
  function(y){
    M = matrix(0,length(y),16)
    sqrty = sqrt(y)
    logy = log(y)
    ry = 1 / y
    rsqrty = 1/sqrty
    if(B[[1]]!=0) {M[,1] = B[[1]] *(ry*ry)             }
    if(B[[2]]!=0) {M[,2] = B[[2]] *(ry)                }
    if(B[[3]]!=0) {M[,3] = B[[3]] *(rsqrty)            }
    if(B[[4]]!=0) {M[,4] = B[[4]] *(logy)              }
    if(B[[5]]!=0) {M[,5] = B[[5]] *(sqrty)             }
    if(B[[6]]!=0) {M[,6] = B[[6]] *y                   }
    if(B[[7]]!=0) {M[,7] = B[[7]] *y*y                 }
    if(B[[8]]!=0) {M[,8] = B[[8]] *y*y*y               }
    if(B[[9]]!=0) {M[,9] = B[[9]] *ry*ry*logy          }
    if(B[[10]]!=0){M[,10]= B[[10]]*ry*logy             }
    if(B[[11]]!=0){M[,11]= B[[11]]*(rsqrty)*logy       }
    if(B[[12]]!=0){M[,12]= B[[12]]*(logy^2)            }
    if(B[[13]]!=0){M[,13]= B[[13]]*sqrty*logy          }
    if(B[[14]]!=0){M[,14]= B[[14]]*y*logy              }
    if(B[[15]]!=0){M[,15]= B[[15]]*y*y*logy            }
    if(B[[16]]!=0){M[,16]= B[[16]]*y*y*y*logy          }
    exp(rowSums(M))
  }
}

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
#'@param rr_specs a list that contains the string variables IM, GENDER, and
#'  FUNCTION, and sixteen double variables B1 through B16
#'
#'@return a function object that describes a base relative risk curve
#'

set_rr <- function(rr_specs) {
  IM <- rr_specs[["IM"]]
  GENDER <- rr_specs[["GENDER"]]
  FUNCTION <- rr_specs[["FUNCTION"]]
  betas <- do.call(paste0, list(rep("B", 16), 1:16))
  BETAS <- as.numeric(rr_specs[betas])
  fn_rr <- function(x) {0}

  if(FUNCTION == "FP"){
    fn_rr <- fp_rr(BETAS)
  }

  if(FUNCTION == "Step" & GENDER == "Female") {
    fn_rr <- hiv_f_rr
  }
  if(FUNCTION == "Step" & GENDER == "Male") {
    fn_rr <- hiv_m_rr
  }

  if(FUNCTION == "Spline") {
    if(IM == "(6).(3)") {
      fn_rr <- acute_pancreatitis_f_rr
    }
    else if(GENDER == "Female") {
      fn_rr <- hypertension_f_rr
    }
    else {
      fn_rr <- hypertension_m_rr
    }
  }
  force(fn_rr)
  fn_rr
}

#' Produce an extrapolation-modified relative risk curve
#'
#'@description given a relative risk curve and some metadata, produces a
#'  risk curve that has the desired extrapolation behaviour beyond 150 g/day,
#'  or 100 g/day for IHD.
#'
#'@param rr_specs a list that contains the string variables IM and GENDER, and
#'  the logical EXT
#'@param rr_fn a relative risk curve
#'
#'@return a function whose values above 150 (100 for IHD) are linearly
#'  extrapolated (see InterMAHP guide for details)
#'

ext_rr <- function(rr_specs, fn_rr) {
  IM <- rr_specs[["IM"]]
  GENDER <- rr_specs[["GENDER"]]
  EXT <- rr_specs[["EXT"]]
  if(GENDER == "Female" & (IM %in% c("(5).(1)", "(6).(3)"))) {
    EXT <- FALSE
  }

  x1 <- 100 + ifelse(IM == "(5).(2)", -50, 0)
  x2 <- 150 + ifelse(IM == "(5).(2)", -50, 0)
  y1 <- fn_rr(x1)
  y2 <- fn_rr(x2)
  slope <- ifelse(EXT, (y2-y1)/(x2-x1), 0)

  lin_ext <- function(x) {
    y2 + slope*(x-x2)
  }

  extrapolated <- function(x) {
    ((0 < x) & (x < x2))*fn_rr(x) + (x2 < x)*lin_ext(x)
  }

  extrapolated
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
#'@param rr_specs a list that contains the string variables IM and the double
#'  BINGEF
#'@param rr_fn a relative risk curve
#'

bng_rr <- function(rr_specs, fn_rr) {
  IM <- rr_specs[["IM"]]
  BINGEF <- rr_specs[["IM"]]
  tmp_rr <- fn_rr
  if(IM %in% c("(5).(2)","(5).(5)")) {
    tmp_rr <- function(x) pmax(1, fn_rr(x))
  }
  bd_rr <- function(x) BINGEF*tmp_rr(x)
  bd_rr
}

#' Compile a list of RR functions
#'
#' @param R a list that contains the string variables IM, GENDER, and FUNCTION,
#' sixteen double variables B1 through B16, the logical EXT, and the double
#' BINGEF.
#'

compile_rr <- function(R) {
  list("BASE_RR" = set_rr(R),
       "LNXT_RR" = ext_rr(R, set_rr(R)),
       "BNGD_RR" = bng_rr(R, set_rr(R)))
}

#' Default Relative Risk Curve Defintion
#'
#' Default input to generate relative risk curves.  Standard formatting is
#' described in the InterMAHP user guides
#'
#'@docType data
#'@usage data(rr_default)
#'
#'
"rr_default"

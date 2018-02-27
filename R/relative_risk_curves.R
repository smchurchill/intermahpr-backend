#' Special Female HIV Relative Risk Function
#'
#'@param x A vector of x values at which to evaluate the RR function
#'
#'
#'
#'

hiv_f_rr <- function(x) {
  ((0<x)&(x<49))*1 + (x>=49)*1.54
}

#' Special Male HIV Relative Risk Function
#'
#'@param x A vector of x values at which to evaluate the RR function
#'
#'
#'
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

#' Get Simple Fractional Polynomial Relative Risk Function
#'
#'@param B The numeric vector of Beta values needed to produce a fractional poly
#'         nomial
#'@param extrapolation Either TRUE(linear) or FALSE(capped) used for
#'                     extrapolating the RR after 150
#'
#'
#'
#'


simple_rr <- function(B, extrapolation = TRUE) {
  force(B)
  force(extrapolation)
  function(x){
    y <- append(x, c(100,150))
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
    sums = exp(rowSums(M))
    x_eval = sums[1:length(x)]
    rr_100 = sums[length(x)+1]
    rr_150 = sums[length(x)+2]
    return(((0 < x) & (x <= 150))*x_eval +
            (x > 150)*(rr_150 +
                       (extrapolation=TRUE)*(x-150)*((rr_150 - rr_100)*0.02)))
  }
}

#' Get Special Fractional Polynomial Relative Risk Function (IHD extrapolated af
#' ter 100)
#'
#'@param B The numeric vector of Beta values needed to produce a fractional poly
#'         nomial
#'@param extrapolation Either TRUE(linear) or FALSE(capped) used for
#'                     extrapolating the RR after 150
#'
#'
#'
#'

ihd_rr <- function(B, extrapolation=TRUE) {
  force(B)
  force(extrapolation)
  function(x){
    y <- append(x, c(50,100))
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
    sums = exp(rowSums(M))
    x_eval = sums[1:length(x)]
    rr_50 = sums[length(x)+1]
    rr_100 = sums[length(x)+2]
    return(((0 < x) & (x <= 100))*x_eval +
            (x>100)*(rr_100 +
                     (extrapolation=TRUE)*(x-100)*(rr_100 - rr_50)*0.02))
  }
}



#' Relative Risk Function Chooser
#'
#'@description Given a row with the standard Relative Risk input structure (plus
#'             extrapolation T/F column), determines a relative risk function
#'             and a relative risk for binge-drinkers function to be returned.
#'
#'@param row is a row from the Relative Risk input frame plus an extrapolation
#'        column.
#'        The necessary structure for each row is
#'        IM > Condition > Gender > Outcome > RR_FD > BingeF > Function > B1-B16
#'        > extrapolation
#'
#'@return A list of:
#'            1. list: IM, Condition Name, Gender, Type (Morb, Mort, Comb)
#'            2. RR_FD constant
#'            3. Relative risk function
#'            4. Relative risk for bingers
#'
#'

rr_chooser <- function(row) {
  function_rr <- function() {0}
  function_bd <- function() {0}

  # Set RR function
  if(row["Function"] == "FP"){
    betas = as.numeric(data.matrix(row[8:23]))
    force(betas)
    if(row[[1]] != "(5).(2)"){
      function_rr <- simple_rr(betas,row[[24]])
    }
    else {
      function_rr <- ihd_rr(betas,row[[24]])
    }

  }
  else if(row[[7]] == "Step") {
    if(row[[3]] == "Female") { function_rr <- hiv_f_rr }
    else { function_rr <- hiv_m_rr }
  }
  else {
    if(row[[1]] == "(6).(3)"){
      function_rr <- acute_pancreatitis_f_rr
    }
    else if(row[[3]] == "Female") {
      function_rr <- hypertension_f_rr
    }
    else {
      function_rr <- hypertension_m_rr
    }
  }

  # Set BD function
  if(row[[6]] != ".") {
    function_bd <- function(x) as.numeric(data.matrix(row[[6]]))*function_rr(x)
  }
  else if(row[[1]] == "(5).(2)" | row[[1]] == "(5).(6)"){
    function_bd <- function(x) pmax(1,function_rr(x))
  }
  else {
    function_bd <- function_rr
  }

  force(function_rr)
  force(function_bd)
  list(row[1:4], row[5], function_rr, function_bd)
}



#' Get a list of Former Drinker Relative Risks, Relative Risk Curves, and Binge
#' Risk Curves
#'
#'@param RR A data table of relative risk information formatted as demanded in
#'           the InterMAHP comprehensive guide
#'
#'@param extrapolation Either TRUE(linear) or FALSE(capped) used for
#'                     extrapolating the RR after 150 (100 for IHD)
#'
#'
#'@export

deduce_relative_risk_curves_from_rr <- function(RR = intermahpr::rr_default,
                                                x_in = TRUE) {
    RR$extrapolation <- x_in
    apply(RR, 1, rr_chooser)
}

#' Default Relative Risk Curve Defintion
#'
#' Default input to generate relative risk curves.  Standard formatting is descr
#' ibed in the InterMAHP user guides
#'
#'@docType data
#'@usage data(rr_default)
#'
#'
"rr_default"

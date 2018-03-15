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

#' Get Fractional Polynomial Relative Risk Function
#'
#'@param B The numeric vector of Beta values needed to produce a fractional poly
#'         nomial
#'
#'
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
                     (extrapolation==TRUE)*(x-100)*(rr_100 - rr_50)*0.02))
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
  IM <- row[["IM"]]
  Gender <- row[["Gender"]]
  Function <- row[["Function"]]
  Outcome <- row[["Outcome"]]
  extrapolation <- row[["extrapolation"]]
  BingeF <- row[["BingeF"]]
  RR_FD <- as.numeric(row[["RR_FD"]])

  tmp_rr <- function() {0}
  betas <- data.matrix(rep(0,16))
  extrapolate_using <- c(100, 150) + ifelse(IM == "(5).(2)", -50, 0)

  # Set RR function
  if(Function == "FP"){
    # print("Function == FP")
    betas <- data.matrix(sapply(row[9:length(row)-1], as.numeric))
    tmp_rr <- fp_rr(betas)
  }

  if(Function == "Step" & Gender == "Female") {
    # print("Function == Step & Gender == Female")
    tmp_rr <- hiv_f_rr
  }
  if(Function == "Step" & Gender == "Male") {
    # print("Function == Step & Gender == Male")
    tmp_rr <- hiv_m_rr
  }

  if(Function == "Spline") {
    # print("Function == Spline")
    if(IM == "(6).(3)") {
      # print("Function == Spline & IM == 6.3")
      tmp_rr <- acute_pancreatitis_f_rr
      extrapolation = FALSE
    }
    else if(Gender == "Female") {
      # print("Function == Spline & IM != 6.3 & Gender == Female")
      tmp_rr <- hypertension_f_rr
      extrapolation = FALSE
    }
    else {
      tmp_rr <- hypertension_m_rr
    }
  }

  fn_rr <- function(x) {
    x1 <- extrapolate_using[[1]]
    x2 <- extrapolate_using[[2]]
    y1 <- tmp_rr(x1)
    y2 <- tmp_rr(x2)
    slope <- ifelse(extrapolation, (y2-y1)/(x2-x1), 0)
    return(((0 < x) & (x < x2))*tmp_rr(x) + (x2 < x)*(y2 + slope*(x-x2)))
  }


  # Set BD function
  fn_bd <- fn_rr
  if(BingeF != ".") {
    fn_bd <- function(x) as.numeric(BingeF)*tmp_rr(x)
  }
  if(IM == "(5).(2)" | IM == "(5).(6)"){
    fn_bd <- function(x) pmax(1,tmp_rr(x))
  }

  # print(ggplot(data = data.frame(x=0), mapping = aes(x=x)) +
  #          stat_function(fun = fn_rr) +
  #          xlim(0.03, 250))
  force(fn_rr)
  force(fn_bd)
  list(row[1:4],RR_FD, fn_rr, fn_bd)
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
#'
#'

deduce_relative_risk_curves_from_rr <- function(RR = intermahpr::rr_default,
                                                x_in = TRUE) {
    RR$extrapolation <- x_in
    apply(RR, 1, rr_chooser)
}

generate_relative_risk_plots <- function(RR = intermahpr::rr_default,
                                         x_in = TRUE) {
  curves <- deduce_relative_risk_curves_from_rr(RR, x_in)
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

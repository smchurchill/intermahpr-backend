#' Combine curve and prevalence/consumption data into an aaf
#'
#'
#'@param C is an entry in the list curve of the form:
#'   [[1]]: a list [IM, Condition, Gender, Outcome]
#'   [[2]]: Relative risk for former drinkers, double
#'   [[3]]: Relative risk for current nonbinging drinkers, function
#'   [[4]]: Relative risk for current binge drinkers, function
#'
#'@param P is a row from the data.table params with entries:
#'   Region  Year    Gender  Age_Group Population
#'   PCC_g/d Cor.Fac Rel.Con P_LA
#'
#' We don't need those first nine to define integrands.  The remaining are:
#'   P_FD        P_CD    P_BD    drinker PCC.drk
#'   Binge bndry G.shape G.scale nc      df
#'   p_bat       R1      R2
#'
#'

compute_aaf <- function(C, P) {
  PX <- sapply(data.frame(t(P[-(1:9)])), function(f) as.numeric(as.character(f)))

  dfgamma <- function(x) {
    PX[["df"]]*dgamma(x, shape=PX[["Gamma_shape"]], scale=PX[["Gamma_scale"]])
  }

  PR <- function(x) {
    dfgamma(x) * (C[[3]](x)-1)
  }
  PB <- function(x) {
    dfgamma(x) * (C[[4]](x)-1)
  }

  fd_numerator <- PX[["P_FD"]] * (as.numeric(C[[2]]) - 1)

  integrand <- function(x) {
    (x          <= PX[["BB"]]) * (PX[["R1"]]*PR(x) + PX[["R2"]]*PB(x)) +
    (PX[["BB"]] <  x         ) * PB(x)
  }

  aaf_dint <- integrate(integrand,
                        intermahpr_constants$lower_bound,
                        intermahpr_constants$upper_bound)$value

  aaf_denominator <- 1 + fd_numerator + aaf_dint

  aaf_numerator <- fd_numerator

  if(PX[["LB"]] < PX[["UB"]]) {
    aaf_numerator <- integrate(integrand, PX[["LB"]], PX[["UB"]])$value
  }

  aaf <- aaf_numerator / aaf_denominator

  aaf_fd <- fd_numerator / aaf_denominator

  list(Region = P[["Region"]],
       Year = P[["Year"]],
       Gender = C[[1]][[3]],
       Age_Group = P[["Age_Group"]],
       IM = C[[1]][[1]],
       Condition = C[[1]][[2]],
       Outcome = C[[1]][[4]],
       FD_Numerator = fd_numerator,
       AAF_Numerator = aaf_numerator,
       AAF_Denominator = aaf_denominator,
       AAF = aaf,
       AAF_FD = aaf_fd)
}

#' Compute AAF numerators over a given interval from a given RR curve / Probabil
#' ity Dist.
#'
#' This method computes, in full generality, the numerator terms found in
#' Section 4.1 of the InterMAHP Comprehensive guide.  Currently unused, as that
#' integral fits on a single line.
#'
#'@param R1 The ratio (P_CD - P_BD) / (P_CD - P_BAT)
#'@param R2 The ratio (P_BD - P_BAT)/ (P_CD - P_BAT)
#'@param PR The function P(x)(RR(X)-1), P is adjusted gamma distribution, RR is
#'          relative risk
#'@param PB The function P(X)(RB(X)-1), P is adjusted gamma distribution, RB is
#'          binge-adjusted relative risk
#'@param bb The binge drinking definition (g/day)
#'@param lb Lower bound of interval
#'@param ub Upper bound of interval
#'
#'

compute_interval_aaf_deprecated <- function(
  R1, R2,
  PR, PB,
  bb, lb, ub
) {
  integrand <- function(x) {
    (x <= bb)*(R1*PR(x) + R2*PB(x)) + (bb < x)*PB(x)
  }
  integrate(integrand, lb, ub)$value
  # integral <- -1
  # if((lb < bb)&&(bb < ub)) {
  #   integral <-
  #     integrate(function(x) ) + , lb, bb)$value
  #     (R1      * integrate(PR,lb,bb)$value) +
  #     (R2      * integrate(PB,lb,bb)$value) +
  #     integrate(PB,bb,ub)$value
  # }
  # else {
  #   integral <-
  #     integrate(function(x) R1*PR(x) + R2*PB(x), lb, ub)$value
  # }
  # return(integral)
}

#' Compute AAF from inputs given as a list L
#'
#'  CUrrently unused
#'
#'


compute_aaf_deprecated <- function(L) {
  integrand <- function(x) {
    (x <= L$bb)*(L$R1*L$PR(x) + L$R2*L$PB(x)) + (L$bb < x)*L$PB(x)
  }
  aaf_cd <- integrate(integrand, L$lb, L$ub)$value
    # compute_interval_aaf(R1 = L$R1,
    #                              R2 = L$R2,
    #                              PR = L$PR,
    #                              PB = L$PB,
    #                              bb = L$bb,
    #                              lb = L$lb,
    #                              ub = L$ub)
  aaf_fd <- L$P_FD * (L$RR_FD - 1)
  aaf_num <- aaf_cd + aaf_fd
  append(
    L[1:7],
    list(AAF = aaf_num / (1 + aaf_num))
  )
}

#' Combine curve and prevalence/consumption data
#'
#'
#'@param C is an entry in the list curve of the form:
#'   [[1]]: a list [IM, Condition, Gender, Outcome]
#'   [[2]]: Relative risk for former drinkers, double
#'   [[3]]: Relative risk for current nonbinging drinkers, function
#'   [[4]]: Relative risk for current binge drinkers, function
#'
#'@param PX is a row from the data.table params with entries:
#'   Region  Year    Gender  Age_Group Population
#'   PCC_l/y Cor.Fac Rel.Con P_LA
#'
#' We don't need those first nine to define integrands.  The remaining are:
#'   P_FD        P_CD    P_BD    drinker PCC.drk
#'   Binge bndry G.shape G.scale nc      df
#'   p_bat       R1      R2



combine <- function(C, PX) {

  P <- as.numeric(PX[-(1:9)])

  PR <- function(x) {
    P[10] * dgamma(x, shape = P[7], scale = P[8]) * (C[[3]](x)-1)
  }
  PB <- function(x) {
    P[10] * dgamma(x, shape = P[7], scale = P[8]) * (C[[4]](x)-1)
  }

  list(
    Region = PX[["Region"]],
    Year = PX[["Year"]],
    Gender = C[[1]][[3]],
    Age_Group = PX[["Age_Group"]],
    IM = C[[1]][[1]],
    Condition = C[[1]][[2]],
    Outcome = C[[1]][[4]],
    R1 = P[12],
    R2 = P[13],
    PR = PR,
    PB = PB,
    bb = P[6],
    P_FD = P[1],
    RR_FD = as.numeric(C[[2]])
    )
}

#' Compute AAF numerators over a given interval from a given RR curve / Probability Dist.
#'
#' This method computes, in full generality, the numerator terms found in
#' Section 4.1 of the InterMAHP Comprehensive guide.
#'
#'@param R1 The ratio (P_CD - P_BD) / (P_CD - P_BAT)
#'@param R2 The ratio (P_BD - P_BAT)/ (P_CD - P_BAT)
#'@param PR The function P(x)(RR(X)-1), P is adjusted gamma distribution, RR is relative risk
#'@param PB The function P(X)(RB(X)-1), P is adjusted gamma distribution, RB is binge-adjusted relative risk
#'@param bb The binge drinking definition (g/day) [[Only used when lb<bb<ub]]
#'@param lb Lower bound of interval
#'@param ub Upper bound of interval
#'
#'

compute_interval_aaf <- function(
  R1, R2,
  PR, PB,
  bb, lb, ub
) {

  integral <- -1
  if((lb < bb)&&(bb < ub)) {
    integral <-
      (R1      * integrate(PR,lb,bb)$value) +
      (R2      * integrate(PB,lb,bb)$value) +
      integrate(PB,bb,ub)$value
  }
  else {
    integral <-
      (R1      * integrate(PR,lb,ub)$value) +
      (R2      * integrate(PB,lb,ub)$value)
  }
  return(integral)
}

#' Compute AAF from inputs given as a list L
#'


compute_aaf <- function(L) {
  aaf_cd <- compute_interval_aaf(R1 = L$R1,
                                 R2 = L$R2,
                                 PR = L$PR,
                                 PB = L$PB,
                                 bb = L$bb,
                                 lb = intermahpr_constants$lower_bound,
                                 ub = intermahpr_constants$upper_bound)
  aaf_fd <- L$P_FD * (L$RR_FD - 1)
  aaf_num <- aaf_cd + aaf_fd
  append(
    L[1:7],
    list(AAF = aaf_num / (1 + aaf_num))
  )
}


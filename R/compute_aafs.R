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
#'
#'@export



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
#'@param bb The binge drinking definition (g/day)
#'@param lb Lower bound of interval
#'@param ub Upper bound of interval
#'
#'@export

compute_interval_aaf <- function(
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
#'@export


compute_aaf <- function(L, lb = intermahpr_constants$lower_bound, ub = intermahpr_constants$upper_bound) {
  aaf_cd <- compute_interval_aaf(R1 = L$R1,
                                 R2 = L$R2,
                                 PR = L$PR,
                                 PB = L$PB,
                                 bb = L$bb,
                                 lb = lb,
                                 ub = ub)
  aaf_fd <- L$P_FD * (L$RR_FD - 1)
  aaf_num <- aaf_cd + aaf_fd
  append(
    L[1:7],
    list(AAF = aaf_num / (1 + aaf_num))
  )
}




#'  Compute Total AAFs
#'
#'@examples
#'  compute_total_aafs(intermahpr::pc_default,
#'                     intermahpr::rr_default,
#'                     TRUE,
#'                     250,
#'                     c(53.80, 67.25))
#'
#'@export
#'
#'
#'@param PC Prevalence and consumption
#'@param RR Relative risk informer
#'@return data.frame of aafs by region, year, gender, age group, condition
#'
#'@description
#'  Compute them.
#'

compute_total_aafs <- function(prev_cons_data = intermahpr::pc_default,
                               relative_risks = intermahpr::rr_default,
                               extrapolation = TRUE,
                               upper_bound = 250,
                               binge_levels = c(53.80, 67.25)) {
  set_upper_bound(upper_bound)
  set_female_binge_barrier(binge_levels[1])
  set_male_binge_barrier(binge_levels[2])

  params <- derive_params_from_pc(prev_cons_data)
  curves <- deduce_relative_risk_curves_from_rr(relative_risks, extrapolation)

  aaf_arguments <- lapply(curves,
                          function(C) apply(params[Gender == C[[1]][[3]], ],
                                            1,
                                            function(P) combine(C,P)))

  aaf_list <- lapply(aaf_arguments, function(x) lapply(x, compute_aaf))

  aaf_frames <- lapply(aaf_list,
                       function(x) data.frame(matrix(unlist(x),
                                                     nrow = length(x),
                                                     byrow = TRUE)))

  aaf_frame <- do.call("rbind", aaf_frames)

  plyr::rename(aaf_frame, c("X1" = "Region",
                            "X2" = "Year",
                            "X3" = "Gender",
                            "X4" = "Age_Group",
                            "X5" = "IM",
                            "X6" = "Condition",
                            "X7" = "Outcome",
                            "X8" = "AAF_Total"))
}


#'  Compute Light AAFs
#'
#'@examples
#'  compute_light_aafs(intermahpr::pc_default, intermahpr::rr_default, TRUE, c())
#'
#'@export
#'
#'
#'@param PC Prevalence and consumption
#'@param RR Relative risk informer
#'@return data.frame of aafs by region, year, gender, age group, condition
#'
#'@description
#'  Compute them.
#'

compute_light_aafs <- function(prev_cons_data = intermahpr::pc_default,
                               relative_risks = intermahpr::rr_default,
                               extrapolation = TRUE,
                               upper_bound = 250,
                               binge_levels = c(53.80, 67.25),
                               lm_barriers = c(10,10)) {
  set_upper_bound(upper_bound)
  set_female_binge_barrier(binge_levels[1])
  set_male_binge_barrier(binge_levels[2])
  set_female_light_moderate_barrier(lm_barriers[1])
  set_male_light_moderate_barrier(lm_barriers[2])

  params <- derive_params_from_pc(PC)
  curves <- deduce_relative_risk_curves_from_rr(RR, extrapolation)

  aaf_arguments <- lapply(curves,
                          function(C) apply(params[Gender == C[[1]][[3]], ],
                                            1,
                                            function(P) combine(C,P)))

  aaf_list <- lapply(aaf_arguments, function(x) lapply(x, compute_aaf))

  aaf_frames <- lapply(aaf_list,
                       function(x) data.frame(matrix(unlist(x),
                                                     nrow = length(x),
                                                     byrow = TRUE)))

  aaf_frame <- do.call("rbind", aaf_frames)

  plyr::rename(aaf_frame, c("X1" = "Region",
                            "X2" = "Year",
                            "X3" = "Gender",
                            "X4" = "Age_Group",
                            "X5" = "IM",
                            "X6" = "Condition",
                            "X7" = "Outcome",
                            "X8" = "AAF_Total"))
}



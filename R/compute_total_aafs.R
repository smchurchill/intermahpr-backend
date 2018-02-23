#'  Compute Total AAFs
#'
#'@examples
#'  compute_total_aafs(intermahpr::pc_default, intermahpr::rr_default)
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

compute_total_aafs <- function(PC = get_prevalence_consumption_data(), RR = get_relative_risk_data()) {
  params <- derive_params_from_pc(PC)
  curves <- deduce_relative_risk_curves_from_rr(RR)

  aaf_arguments <- lapply(curves,
                          function(C) apply(params[Gender == C[[1]][[3]], ],
                                            1,
                                            function(P) combine(C,P)))

  aaf_list <- lapply(aaf_arguments, function(x) lapply(x, compute_aaf))

  aaf_frames <- lapply(aaf_list,
                       function(x) data.frame(matrix(unlist(x),
                                                     nrow = 6,
                                                     byrow = T)))

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

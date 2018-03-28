#' Collect and assemble AAF data from RR and PC input data
#'
#'@param pc  Prevalence / Consumption input
#'@param rr  Relative Risk input
#'@param ext logical, extrapolate linearly?
#'@param lb  Double, consumption lower bound
#'@param ub  Double, consumption upper bound
#'@param bb  Double vector, Binge consumption level, Gender stratified
#'
#'

assemble <- function(pc, rr, ext, lb, ub, bb,
                     gc = list("Female" = 1.258^2, "Male" = 1.171^2)) {
  RRF <- format_v0_rr(rr = rr)
  PCF <- format_v0_pc(pc = pc)

  RRD <- derive_v0_rr(rr = RRF, ext = ext)
  PCD <- derive_v0_pc(pc = PCF, bb = bb, lb = lb, ub = ub, gc = gc)

  join_pc_rr(pc = PCD, rr = RRD)
}

#' Compute AAFs for all conditions as separated by the given cutpoints
#'
#'@param aaf_table A tibble as returned by assemble
#'@param cuts A list of double vectors indexed by aaf_table$GENDER
#'
#'

compute_aafs <- function(aaf_table, cuts) {
  aaf_table %<>%
    mutate(CUTS = pmap(list(LB, cuts[GENDER], UB), c)) %>%
    mutate(CUMUL_F = map2(AAF_CMP, CUTS, ~.x(.y))) %>%
    mutate(AAF_GRP = map(CUMUL_F, ~diff(.x))) %>%
    mutate(AAF_TOT = map(CUMUL_F, ~`[[`(.x, length(.x))))
  aaf_table
}




#'  Compute General AAFs
#'
#'@examples
#'  compute_general_aafs(intermahpr::pc_default,
#'                     intermahpr::rr_default,
#'                     TRUE,
#'                     250,
#'                     c(53.80, 67.25))
#'
#'
#'@export
#'
#'@param PC Prevalence and consumption
#'@param RR Relative risk informer
#'@return data.frame of aafs by region, year, gender, age group, condition
#'
#'@description
#'  Interface functions are wrappers for this more general aaf compiler
#'

compute_general_aafs <- function(prev_cons_data = intermahpr::pc_default,
                                 relative_risks = intermahpr::rr_default,
                                 extrapolation = TRUE,
                                 upper_bound = 250,
                                 binge_levels = c(53.80, 67.25),
                                 lower = c(0.03,0.03),
                                 upper = c(upper_bound, upper_bound)) {
  set_upper_bound(upper_bound)

  params <- derive_params_from_pc(prev_cons_data,
                                  bb = binge_levels,
                                  lb = lower,
                                  ub = upper)
  curves <- deduce_relative_risk_curves_from_rr(relative_risks, extrapolation)


  aaf_list <- lapply(curves,
                     function(C) apply(subset(params, Gender == C[[1]][[3]]),
                                       1,
                                       function(x) compute_aaf(C,x)))

  saved_names <- names(aaf_list[[1]][[1]])

  aaf_frames <- lapply(aaf_list,
                       function(x) data.frame(matrix(unlist(x),
                                                     nrow = length(x),
                                                     byrow = TRUE)))


  aaf_frame <- do.call("rbind", aaf_frames)

  names(aaf_frame) <- saved_names

  aaf_frame
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
#'  Compute total aafs
#'

compute_total_aafs <- function(prev_cons_data = intermahpr::pc_default,
                               relative_risks = intermahpr::rr_default,
                               extrapolation = TRUE,
                               upper_bound = 250,
                               binge_levels = c(53.80, 67.25)) {

  aafs <- compute_general_aafs(prev_cons_data, relative_risks, extrapolation,
                               upper_bound, binge_levels,
                               lower = rep(intermahpr_constants$lower_bound, 2),
                               upper = rep(upper_bound, 2))

  aaf_rename <- plyr::rename(aafs[-(8:10)], c("AAF" = "AAF_CD"))

  DT <- data.table(sapply(aaf_rename[-(1:7)],
                          function(f) as.numeric(as.character(f))))
  DT[, AAF_Total := AAF_FD + AAF_CD]
  cbind(aafs[1:7], DT)
}

#'  Compute Light AAFs
#'
#'@examples
#' compute_light_aafs(intermahpr::pc_default, intermahpr::rr_default, TRUE, c())
#'
#'@export
#'
#'
#'@param PC Prevalence and consumption
#'@param RR Relative risk informer
#'@return data.frame of aafs by region, year, gender, age group, condition
#'
#'@description
#'  Compute light aafs
#'

compute_light_aafs <- function(prev_cons_data = intermahpr::pc_default,
                               relative_risks = intermahpr::rr_default,
                               extrapolation = TRUE,
                               upper_bound = 250,
                               binge_levels = c(53.80, 67.25),
                               light_moderate_barriers = c(20,20)) {

  aafs <- compute_general_aafs(prev_cons_data, relative_risks, extrapolation,
                               upper_bound, binge_levels,
                               lower = rep(intermahpr_constants$lower_bound, 2),
                               upper = light_moderate_barriers)

  plyr::rename(aafs[-(8:10)], c("AAF" = "AAF_LD"))

}

#'  Compute Moderate AAFs
#'
#'@examples
#' compute_mod_aafs(intermahpr::pc_default, intermahpr::rr_default, TRUE, c())
#'
#'@export
#'
#'
#'@param PC Prevalence and consumption
#'@param RR Relative risk informer
#'@return data.frame of aafs by region, year, gender, age group, condition
#'
#'@description
#'  Compute moderate aafs
#'

compute_mod_aafs <- function(prev_cons_data = intermahpr::pc_default,
                               relative_risks = intermahpr::rr_default,
                               extrapolation = TRUE,
                               upper_bound = 250,
                               binge_levels = c(53.80, 67.25),
                               light_moderate_barriers = c(20,20),
                               moderate_heavy_barriers = c(40,40)) {

  aafs <- compute_general_aafs(prev_cons_data, relative_risks, extrapolation,
                               upper_bound, binge_levels,
                               lower = light_moderate_barriers,
                               upper = moderate_heavy_barriers)

  plyr::rename(aafs[-(8:10)], c("AAF" = "AAF_MD"))


}

#'  Compute Heavy AAFs
#'
#'@examples
#' compute_heavy_aafs(intermahpr::pc_default, intermahpr::rr_default, TRUE, c())
#'
#'@export
#'
#'
#'@param PC Prevalence and consumption
#'@param RR Relative risk informer
#'@return data.frame of aafs by region, year, gender, age group, condition
#'
#'@description
#'  Compute heavy aafs.
#'

compute_heavy_aafs <- function(prev_cons_data = intermahpr::pc_default,
                               relative_risks = intermahpr::rr_default,
                               extrapolation = TRUE,
                               upper_bound = 250,
                               binge_levels = c(53.80, 67.25),
                               moderate_heavy_barriers = c(40,40)) {


  aafs <- compute_general_aafs(prev_cons_data, relative_risks, extrapolation,
                               upper_bound, binge_levels,
                               lower = moderate_heavy_barriers,
                               upper = rep(upper_bound, 2))

  plyr::rename(aafs[-(8:10)], c("AAF" = "AAF_HD"))
}

#'  Compute All AAFs
#'
#'@examples
#' compute_all_aafs(intermahpr::pc_default, intermahpr::rr_default, TRUE, c())
#'
#'@export
#'
#'
#'@param PC Prevalence and consumption
#'@param RR Relative risk informer
#'@return data.frame of aafs by region, year, gender, age group, condition
#'
#'@description
#'  Compute all AAFs.
#'

compute_all_aafs <- function(prev_cons_data = intermahpr::pc_default,
                             relative_risks = intermahpr::rr_default,
                             extrapolation = TRUE,
                             upper_bound = 250,
                             binge_levels = c(53.80, 67.25),
                             light_moderate_barriers = c(20,20),
                             moderate_heavy_barriers = c(40,40)) {

  # t <- compute_total_aafs(prev_cons_data, relative_risks, extrapolation,
  #                         upper_bound, binge_levels)

  l <- compute_light_aafs(prev_cons_data, relative_risks, extrapolation,
                          upper_bound, binge_levels, light_moderate_barriers)
  m <- compute_mod_aafs(prev_cons_data, relative_risks, extrapolation,
                        upper_bound, binge_levels,
                        light_moderate_barriers, moderate_heavy_barriers)
  h <- compute_heavy_aafs(prev_cons_data, relative_risks, extrapolation,
                          upper_bound, binge_levels, moderate_heavy_barriers)


  all_aafs <- Reduce(function(...) merge(..., all=TRUE), list(l, m, h))
  # , t))
  DT <- data.table(sapply(all_aafs[-(1:7)],
                          function(f) as.numeric(as.character(f))))
  DT[, AAF_Total := AAF_FD + AAF_LD + AAF_MD + AAF_HD]
  # DT[, sanity_checker := computed_total - AAF_Total]

  cbind(all_aafs[1:7], DT)
}




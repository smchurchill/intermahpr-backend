makeCurrentFreeFactory <- function(extrapolated_risk, binge_risk, rr_fd) {
  ## *d_cmp_fct functions have positional ... arguments. positions are:
  ## ..1: population exposure mass function
  ## ..2: proportion of drinkers below BB that do not binge
  ## ..3: proportion of drinkers below BB that do binge
  ## ..4: lower bound
  ## ..5: binge barrier
  ## ..6: upper bound
  ##
  ## for independent current drinker aafs, all are used.
  function(...) {
    former_comp_factory <- makeFormerFreeComponentFactory(rr_fd = rr_fd)
    former_comp <- former_component_factory()
    current_comp_factory <- makeCurrentFreeComponentFactory(
      extrapolated_risk = extrapolated_risk, binge_risk = binge_risk)
    current_comp <- current_component_factory(
      exposure_mass = ..1, non_bingers = ..2, bingers = ..3, lb = ..4, bb = ..5)
    reciprocal_denom <- 1 / (1 + former_comp() + current_comp(..6))
    function(x) reciprocal_denom * current_comp(x)
  }
}

makeCurrentFreeComponentFactory <- function(extrapolated_risk, binge_risk) {
  function(exposure_mass, non_bingers, bingers, lb, bb) {
    pf <- makePreventableFraction(
      bb = bb, r1 = r1, r2 = r2,
      lnxt_rr = lnxt_rr, bngd_rr = bngd_rr
    )
    cdf_fct_fct(f = n_gamma %prod% pf, lb = lb)
  }
}

makeFormerFreeFactory <- function(extrapolated_risk, binge_risk, rr_fd) {
  ## *d_cmp_fct functions have positional ... arguments. positions are:
  ## ..1: population exposure mass function
  ## ..2: proportion of drinkers below BB that do not binge
  ## ..3: proportion of drinkers below BB that do binge
  ## ..4: lower bound
  ## ..5: binge barrier
  ## ..6: upper bound
  ##
  ## for independent former drinker aafs, all are used.
  function(...) {
    fd_comp_f <- makeFormerFreeComponentFactory(rr_fd = rr_fd)
    fd_comp <- fd_comp()
    cd_comp_f <- makeCurrentFreeComponentFactory(lnxt_rr = lnxt_rr, bngd_rr = bngd_rr)
    cd_comp <- cd_comp_f(n_gamma = ..1, r1 = ..2, r2 = ..3, lb = ..4, bb = ..5)
    rdenom <- 1 / (1 + fd_comp() + cd_comp(..6))
    function(x) rdenom * fd_comp
  }
}

makeFormerFreeComponentFactory <- function(rr_fd){
  function(p_fd) {
    function(x) {
      p_fd*(rr_fr-1)
    }
  }
}

#' Factory for preventable fractions
#'
#'@description Produces the combined and scaled function that respresents the
#'  preventable fraction of disease that, when integrated against exposure,
#'  produces an attributable fraction term.
#'
#'
#'@param bb dbl, binge barrier
#'@param r1 dbl, proportion of drinkers below BB that do not binge
#'@param r2 dbl, proportion of drinkers below BB that do binge
#'@param lnxt_rr fn, extrapolated relative risk for nonbingers
#'@param bngd_rr fn, extrapolated relative risk for bingers
#'
#'

makePreventableFraction <- function(bb, non_binge, binge, ext_risk, binge_risk){
  function(x) {
    (x<=bb)*(r1*(lnxt_rr(x)-1) + r2*(bngd_rr(x)-1)) + (x>bb)*(bngd_rr(x)-1)
  }
}

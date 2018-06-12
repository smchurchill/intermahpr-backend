#### Well-defined AAF Computer Factories ---------------------------------------

#' Factory for current drinker's AAF computer factory for a condition with
#'  well-defined relative risk
#'

makeCurrentFreeFactory <- function(ext_risk, binge_risk, rr_fd) {
  force(ext_risk)
  function(args) {
    former_comp_factory <- makeFormerFreeComponentFactory(rr_fd = rr_fd)
    former_comp <- former_comp_factory(p_fd = args$p_fd)
    current_comp_factory <- makeCurrentFreeComponentFactory(
      ext_risk = ext_risk,
      binge_risk = binge_risk)
    current_comp <- current_comp_factory(args)
    reciprocal_denom <- 1 / (1 + former_comp() + current_comp(args$ub))
    function(x) reciprocal_denom * current_comp(x)
  }
}

#' Factory for former drinker's AAF computer factory for a condition with
#' well-defined relative risk
#'

makeFormerFreeFactory <- function(ext_risk, binge_risk, rr_fd) {
  function(args) {
    former_comp_factory <- makeFormerFreeComponentFactory(rr_fd = rr_fd)
    former_comp <- former_comp_factory(p_fd = args$p_fd)
    current_comp_factory <- makeCurrentFreeComponentFactory(
      ext_risk = ext_risk,
      binge_risk = binge_risk)
    current_comp <- current_comp_factory(args)
    reciprocal_denom <- 1 / (1 + former_comp() + current_comp(args$ub))
    function(x) reciprocal_denom * former_comp()
  }
}

#### AAF Component Factories ---------------------------------------------------

#' Factory for the current drinker's component in an AAF computer factory for a
#'  condition with well-defined relative risk
#'

makeCurrentFreeComponentFactory <- function(ext_risk, binge_risk) {
  function(args) {
    preventable_fraction <- makePreventableFraction(
      bb = args$bb,
      non_bingers = args$non_bingers,
      bingers = args$bingers,
      ext_risk = ext_risk,
      binge_risk = binge_risk)
    integrand <- args$mass %prod% preventable_fraction
    makeIntegrator(f = integrand, lb = args$lb, ub = args$ub)
  }
}

#' Factory for the former drinker's component in an AAF computer factory for a
#'  condition with well-defined relative risk
#'

makeFormerFreeComponentFactory <- function(rr_fd){
  function(p_fd) {
    function(...) {
      p_fd*(rr_fd-1)
    }
  }
}

#' Factory for preventable fraction functions
#'
#'@description Produces the combined and scaled function that respresents the
#'  preventable fraction of disease that, when integrated against exposure,
#'  produces an attributable fraction term.
#'
#'
#'@param bb dbl, binge barrier
#'@param non_bingers dbl, proportion of drinkers below BB that do not binge
#'@param bingers dbl, proportion of drinkers below BB that do binge
#'@param ext_risk fn, extrapolated relative risk for nonbingers
#'@param binge_risk fn, extrapolated relative risk for bingers
#'
#'

makePreventableFraction <- function(
  bb, non_bingers, bingers, ext_risk, binge_risk
) {
  function(x) {
    (x<=bb)*(non_bingers*(ext_risk(x)-1) + bingers*(binge_risk(x)-1)) +
    (x>bb)*(binge_risk(x)-1)
  }
}

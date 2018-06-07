#### Gamma Function Factories --------------------------------------------------

#' Factory for gamma distributions
#'
#'@description
#'We want easier access to gamma distributions with a wide array of fixed params
#'
#'@param sh,sc Gamma parameters to be supplied to dgamma
#'
#'

makeGamma <- function(sh, sc) {
  function(x) dgamma(x, shape = sh, scale = sc)
}

#' Factory for normalized gamma distributions
#'
#'@param shape,scale Input gamma distribution params
#'@param factor Factor to scale GAMMA by
#'
#'

makeNormalizedGamma <- function(shape, scale, factor) {
  function(x) factor * dgamma(x, shape = shape, scale = scale)
}

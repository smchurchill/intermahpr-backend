#' Compute AAFs for all conditions for current drinkers and totals AAFs
#'
#'@param aaf_table A tibble as returned by assemble (i.e. has the following vars
#' AAF_CMP: fn, VIVO, computes AAF from LB to each input
#' AAF_FD: dbl, AAF for former drinkers
#' UB: dbl, consumption upper bound
#')
#'
#'@importFrom purrr pmap map map2
#'@importFrom magrittr "%>%" "%<>%"
#'
#'@export
#'

base_aafs <- function(aaf_table) {
  aaf_table %<>%
    mutate(AAF_CD = unlist(map2(AAF_CMP, UB, ~(.x(.y))))) %>%
    mutate(AAF_TOTAL = AAF_FD + AAF_CD)
  aaf_table
}


#' Compute AAFs for all conditions as separated by the given cutpoints
#'
#'@param aaf_table A tibble as returned by assemble (i.e. has the following vars
#' GENDER: chr, gender levels
#' AAF_CMP: fn, VIVO, computes AAF from LB to each input
#' AAF_FD: dbl, AAF for former drinkers
#' LB: dbl, consumption lower bound
#' UB: dbl, consumption upper bound
#')
#'@param cuts A sorted list of double vectors indexed by aaf_table$GENDER st.
#' each value is between LB and UB
#'
#'@importFrom purrr pmap map map2
#'@importFrom magrittr "%>%" "%<>%"
#'
#'@export
#'

compute_aafs <- function(aaf_table, cuts) {
  aaf_table %<>%
    mutate(CUTS = pmap(list(LB, cuts[GENDER], UB), c)) %>%
    mutate(CUMUL_F = map2(AAF_CMP, CUTS, ~.x(.y))) %>%
    mutate(AAF_GRP = map(CUMUL_F, ~diff(.x))) %>%
    mutate(AAF_CD = unlist(map(CUMUL_F, ~`[`(.x, length(.x))))) %>%
    mutate(AAF_TOTAL = AAF_FD + AAF_CD)
  aaf_table
}

#' Add evaluation cutpoints to a given datatable
#'
#' Written to provide additional reactivity to Shiny InterMAHP app.
#'
#'@param aaf_table a tibble as returned by assemble
#'@param cuts a list of double vectors indexed by aaf_table$GENDER
#'
#'@importFrom purrr pmap
#'@importFrom magrittr "%<>%"
#'
#'@export
#'

add_cutpoints <- function(aaf_table, cuts) {
  aaf_table %<>%
    mutate(CUTS = pmap(list(LB, cuts[GENDER], UB), c))
  aaf_table
}

#' Evaluate a given datatable at pre-added cutpoints
#'
#' Written to provide additional reactivity to shiny app
#'
#'@param aaf_table_cuts a tibble as returned by add_cutpoints
#'
#'@importFrom purrr pmap map map2
#'@importFrom magrittr "%>%" "%<>%"
#'
#'@export
#'

evaluate_at_cutpoints <- function(aaf_table_cuts) {
  aaf_table_cuts %<>%
    mutate(CUMUL_F = map2(AAF_CMP, CUTS, ~.x(.y))) %>%
    mutate(AAF_GRP = map(CUMUL_F, ~diff(.x))) %>%
    mutate(AAF_CD = unlist(map(CUMUL_F, ~`[`(.x, length(.x))))) %>%
    mutate(AAF_TOTAL = AAF_FD + AAF_CD)
  aaf_table_cuts
}

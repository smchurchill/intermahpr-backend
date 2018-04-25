#' Compute AAFs for all conditions as separated by the given cutpoints
#'
#'@param aaf_table A tibble as returned by assemble
#'@param cuts A list of double vectors indexed by aaf_table$GENDER
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

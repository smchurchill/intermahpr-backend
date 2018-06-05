#### Process DH and PC data ----------------------------------------------------

#' Joins DH and PC data
#'
#'@param dh Deaths/Hosps as returned by format_dh
#'@param pc Prev/Cons as returned by derive_pc
#'
#'@importFrom magrittr %>% %<>%
#'
#'@export

process_dh_pc <- function(dh, pc) {
  dplyr::inner_join(
    dh,
    pc[c("REGION", "YEAR", "GENDER", "AGE_GROUP",
         "DRINKERS", "BB", "LB", "UB", "N_GAMMA")],
    by = c("REGION", "YEAR", "GENDER", "AGE_GROUP")
  )
}


#### Process AAF and DH Data ---------------------------------------------------

#' Joins DH and AAF data, fills in missing aaf_cmps
#'
#'@description
#'Provides continuous AAF computing functions for 100% attributable conditions
#'either by calibrating from counts and drinkers (3.2|4.[123]|5.3|6.[15]) or
#'rescaling similar aaf_cmps (4.[467]|8.5|9.2).
#'Granulates [89].all into respective boxes, copies 6.2 into 5.7.
#'
#'@param dh death/hosp tibble as produced by derive_dh
#'@param aaf aaf tibble as produced by base_aafs
#'
#'
#'@importFrom magrittr %>% %<>%
#'
#'@export
#'


process_dh_aaf <- function(dh, aaf) {
  ## Variables to join DH and AAF by
  SIMILAR <- c("IM", "REGION", "YEAR", "GENDER", "AGE_GROUP", "OUTCOME")

  ## [89].all will be distributed and dropped, so don't need as many vars
  NO_INJ <- aaf %>%
    filter(!grepl("all", IM)) %>%
    select(-c(CONDITION, DRINKERS, BB, LB, UB, N_GAMMA))

  INJ <- aaf %>%
    filter(grepl("all", IM)) %>%
    mutate(COUNT = -1)

  ## Join the tables and add Block and Category variables
  JOIN <- dplyr::bind_rows(
    left_join(x = dh, y = NO_INJ, by = SIMILAR),
    INJ
    ) %>%
    mutate(BLOCK = as_factor(purrr::pmap_chr(
      list(REGION, YEAR, GENDER, AGE_GROUP, OUTCOME),
      .f = paste0))
    ) %>%
    mutate(CC = substr(IM, 1, 3))

  ## These rows are fine as is
  KEEP <- filter(JOIN,!grepl("(3...2|4...[^5]|5...[37]|6...[15]|^.[89])", IM))

  ## These rows get a calibrated AAF_CMP
  CLBR <- JOIN %>%
    filter(
      grepl("(3...2|4.*[123]|5.*3|6.*[15])", IM)
    ) %>%
    mutate(
      AAF_CMP = pmap(
        list(IM, COUNT, DRINKERS, N_GAMMA, LB, BB, UB),
        aaf_calibration_factory
      ),
      AAF_FD = 0
    )

  ## These rows get rescaled AAF_CMPs from similar conditions
  SIMS_L <- JOIN %>%
    filter(
      grepl("(4.*[467]|8...5|9...2)", IM)
    ) %>%
    select(-c(AAF_FD, AAF_CMP, AAF_TOTAL))
  SIMS_R <- JOIN %>%
    filter(
      grepl("(4.*5|all)", IM)
    ) %>%
    select(c(AAF_FD, AAF_CMP, AAF_TOTAL, BLOCK, CC))

  SIMS <- left_join(x = SIMS_L, y = SIMS_R, by = c("BLOCK", "CC")) %>%
    mutate(
      AAF_CMP = map2(1 / AAF_TOTAL, AAF_CMP, ~{function(x) .x * .y(x)}),
      AAF_TOTAL = 1)


  ## These rows get full copies of AAF_CMP from related conditions
  COPY_L <- JOIN %>%
    filter(
      grepl("(5.*7|8...[12346]|9...[1345])", IM)
    ) %>%
    select(-c(AAF_FD, AAF_CMP, AAF_TOTAL))
  COPY_R <- JOIN %>%
    filter(
      grepl("(6.*2|all)", IM)
    ) %>%
    select(
      c(AAF_FD, AAF_CMP, AAF_TOTAL, BLOCK, CC)
    ) %>%
    mutate(
      CC = ifelse(CC == "(6)", "(5)", CC)
    )

  COPY <- left_join(x = COPY_L, y = COPY_R, by = c("BLOCK", "CC"))

  ## Combine everything back together, recompute total aafs, and apply to counts
  bind_rows(KEEP, CLBR, SIMS, COPY) %>%
    mutate(AAF_CD = map2_dbl(AAF_CMP, UB, ~(.x(.y)))) %>%
    mutate(AAF_TOTAL = AAF_CD + AAF_FD) %>%
    mutate(AA_CD_COUNT = AAF_TOTAL * COUNT) %>%
    mutate(AA_FD_COUNT = AAF_FD * COUNT)
}

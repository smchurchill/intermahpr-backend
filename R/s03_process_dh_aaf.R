#### Process AAF and DH Data ---------------------------------------------------

#' Joins DH and AAF data, fills in missing aaf_cmps
#'
#'@description
#'Provides continuous AAF computing functions for 100% attributable conditions
#'either by calibrating from counts and drinkers (3.2|4.[123]|5.3|6.[15]) or
#'rescaling similar aaf_cmps (4.[467]|8.5|9.2).
#'Granulates [89].all into respective boxes, copies 6.2 into 5.7.
#'
#'@param dh death/hosp tibble as produced by derive_v*_dh
#'@param aaf aaf tibble as produced by base_aafs
#'
#'
#'@importFrom magrittr %>% %<>%
#'
#'@export
#'


join_dh_aaf <- function(dh, aaf) {
  SIMILAR <- c("IM", "REGION", "YEAR", "GENDER", "AGE_GROUP", "OUTCOME")

  NO_INJ <- aaf %>%
    filter(!grepl("all", IM)) %>%
    select(-c(CONDITION, DRINKERS, BB, LB, UB, N_GAMMA))

  INJ <- aaf %>%
    filter(grepl("all", IM)) %>%
    mutate(COUNT = -1)

  JOINED <- dplyr::bind_rows(
    left_join(x = dh, y = NO_INJ, by = SIMILAR),
    INJ
    ) %>%
    mutate(BLOCK = as_factor(purrr::pmap_chr(
      list(REGION, YEAR, GENDER, AGE_GROUP, OUTCOME),
      .f = paste0))
    ) %>%
    mutate(KEY = as_factor(purrr::pmap_chr(
      list(IM, REGION, YEAR, GENDER, AGE_GROUP, OUTCOME),
      .f = paste0))
    ) %>%
    mutate(CC = substr(IM, 1, 3))

  KEEP <- filter(JOINED,!grepl("(4...[^5]|5...[37]|6...[15]|[89]..\\()", IM))
  MOD  <- filter(JOINED, grepl("(4...[^5]|5...[37]|6...[15]|[89]..\\([^a])", IM))
  USE  <- filter(JOINED, grepl("(4.*5|6.*2|all)", IM))

  for(n in 1:nrow(MOD)) {
    im <- MOD[[n, "IM"]]
    if(grepl("(4...[123]|5...3|6...[15])", im)) {
      expr <- NULL
      use <- MOD[n, ]
      get_aaf_fn <- calibration_factory
    }
    else {
      block <- MOD[[n, "BLOCK"]]
      cc <- MOD[[n, "CC"]]
      if(grepl("(4...[467]|8...5|9...2)", im)) {
        expr <- "(4.*5|all)"
        get_aaf_fn <- scaled_aaf_cmp_factory
      }
      else if(grepl("(5...7|8...[^5]|9...[^2])", im)) {
        if(grepl("5...7", im)) {
          cc <- "(6)"
        }
        expr <- "(6.*2|all)"
        get_aaf_fn <- copy_aaf_cmp_factory
      }
      else {
        stop(
          paste0(
            "IM of ",
            im,
            " does not match any of the following regular expressions:\n",
            "(4...[123]|5...3|6...[15])\n",
            "(4...[467]|8...5|9...2)\n",
            "(5...7|8...[^5]|9...[^2])",
            collapse = "\n"
          )
        )
      }
      use <- filter(USE, grepl(expr, IM) & BLOCK == block & CC == cc)
    }
    MOD[[n, "AAF_CMP"]] <- get_aaf_fn(use)
    MOD[[n, "AAF_FD"]] <- aaf_fd(use)
  }

  bind_rows(KEEP, MOD)
}

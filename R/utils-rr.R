dep_rr <- function(rr, pc, dh) {
  .data <- inner_join(
    x = rr,
    y = inner_join(
      x = pc, y = dh, by = c(region, year, gender, age_group, outcome)
    ),
    by = c(im, gender, outcome)
  ) %>%
    mutate(
      threshold = ifelse(grepl("6", im), lb, bb),
      incidence = count / drinkers
    ) %>%
    mutate(
      aaf_cmp_fct = pmap(
        list(incidence, n_gamma, threshold, ub),
        dep_aaf_cmp_fct_fct
      )
    )
}


ind_rr <- function() {}

filter_dep <- function(.data) {
  filter(.data, outcome == "Calibrated")
}

filter_ind <- function(.data) {
  filter(.data, outcome != "Calibrated")
}


#' Crushes numbered variables
#'
#' All numbered variables (i.e. b1-b16) are crushed into a single variable.
#' Assumes numbered variables are al of numeric type.
#' New variable is named 'betas'.
#' variables with internal numbers (i.e. free2delete) are removed but not
#' crushed
#'

crush_betas <- function(.data) {
  .crush <- .data[grep("[0-9]$", names(.data))]
  .data <- .data[-grep("[0-9]", names(.data))]

  .crushed <- split(as.matrix(.crush), 1:nrow(.crush))
  .data$betas <- .crushed
  .data
}


#' List of variables expected to be in an RR sheet
#'
#'

rr_vars <- c(
  "im",
  "condition",
  "gender",
  "outcome",
  "rr_fd",
  "bingef",
  "form",
  "attributability",
  paste0("b", 1:16)
)

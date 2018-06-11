#### Relative Risk Specific Data Carpentry -------------------------------------

#' Prepare Relative Risk Data
#'

prepareRR <- function(.data, ext) {
  .data %<>%
    cleanRR() %>%
    mutate(ext = ext) %>%
    splitOutcome() %>%
    splitGender()
}

#' Clean Relative Risk Data

cleanRR <- function(.data) {
  crushBetas(clean(.data, rr_vars))
}

#' Split 'Combined' outcomes into Morbidity and Mortality

splitOutcome <- function(.data) {
  morb <- .data %>%
    filter(grepl("(b|C)", outcome)) %>%
    mutate(outcome = "Morbidity")

  mort <- .data %>%
    filter(grepl("(t|C)", outcome)) %>%
    mutate(outcome = "Mortality")

  rbind(morb, mort)
}

#' Split 'All' genders into 'Female' and 'Male'
#'

splitGender <- function(.data) {
  female <- .data %>%
    filter(grepl("(F|A)", gender)) %>%
    mutate(gender = "Female")

  male <- .data %>%
    filter(grepl("(M|A)", gender)) %>%
    mutate(gender = "Male")

  rbind(female, male)
}

#' Filter calibrated forms from the rest

filterCalibrated <- function(.data) {
  filter(.data, form == "Calibrated")
}

#' Filter well-defined forms from the rest

filterFree <- function(.data) {
  filter(.data, form != "Calibrated")
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

#' Crush numbered variables
#'
#' All numbered variables (i.e. b1-b16) are crushed into a single variable.
#' Assumes numbered variables are al of numeric type.
#' New variable is named 'betas'.
#'

crushBetas <- function(.data) {
  crush <- .data[grep("[0-9]$", names(.data))]
  .data <- .data[-grep("[0-9]", names(.data))]

  crushed <- split(as.matrix(crush), 1:nrow(crush))
  .data$betas <- crushed
  .data
}

#### Factories -----------------------------------------------------------------

#' Factory for AAF computer factories: conditions requiring calibration against
#' population statistics

makeCalibratedFactories <- function(rr, pc, dh) {
  .data <- inner_join(
    x = rr,
    y = inner_join(
      x = pc, y = dh, by = c("region", "year", "gender", "age_group")
    ),
    by = c("im", "gender", "outcome")
    ) %>%
    mutate(
      threshold = ifelse(grepl("6", im), lb, bb),
      incidence = count / drinkers
    ) %>%
    mutate(
      current_fraction_factory = pmap(
        list(target = incidence, clbr_mass = n_gamma, lb = threshold, ub = ub),
        makeCurrentCalibratedFactory
      ),
      former_fraction_factory = pmap(
        list(NA),
        makeFormerCalibratedFactory
      )
    )
}



#' Factory for AAF computer factories: conditions with well-defined rel. risk

makeFreeFactories <- function(.data) {
  .data %>%
    mutate(
      base_risk = pmap(
        list(im  = im, gender = gender, form = form, betas = betas),
        makeBaseRisk
      )
    ) %>%
    mutate(x1 = ifelse(grepl("5...2", im), 50, 100)) %>%
    mutate(x2 = x1 + 50) %>%
    mutate(y1 = map2_dbl(x1, base_risk, ~.y(.x)),
           y2 = map2_dbl(x2, base_risk, ~.y(.x))) %>%
    mutate(slope = ifelse(ext, (y2-y1)/(x2-x1), 0)) %>%
    mutate(
      ext_risk = pmap(
        list(base_risk = base_risk, x2 = x2, y2 = y2, slope = slope),
        makeExtrapolatedRisk
      )
    ) %>%
    mutate(
      binge_risk = pmap(
        list(im = im, bingef = bingef, ext_risk = ext_risk),
        makeBingeRisk
      )
    ) %>%
    mutate(
      current_fraction_factory = pmap(
        list(
          ext_risk = ext_risk,
          binge_risk = binge_risk,
          rr_fd = rr_fd),
        makeCurrentFreeFactory
      ),
      former_fraction_factory = pmap(
        list(
          ext_risk = ext_risk,
          binge_risk = binge_risk,
          rr_fd = rr_fd),
        makeFormerFreeFactory
      )
    )
}






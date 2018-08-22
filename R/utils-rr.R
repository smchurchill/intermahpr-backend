#### Relative Risk Specific Data Carpentry -------------------------------------

#' Prepare Relative Risk Data
#'
#' @export
prepareRR <- function(.data, ext) {
  message("Preparing relative risk input... ", appendLF = FALSE)
  .data %<>%
    clean(getExpectedVars("rr")) %>%
    crushBetas() %>%
    mutate(ext = ext) %>%
    splitOutcome() %>%
    splitGender()

  message("Done")

  .data
}

#' Split 'Combined' and 'Calibrated' outcomes into Morbidity and Mortality
#' @export
splitOutcome <- function(.data) {
  morb <- .data %>%
    filter(grepl("(Morbidity|Calibrated|Combined)", outcome)) %>%
    mutate(outcome = "Morbidity")

  mort <- .data %>%
    filter(grepl("(Mortality|Calibrated|Combined)", outcome)) %>%
    mutate(outcome = "Mortality")

  rbind(morb, mort)
}

#' Split 'All' genders into 'Female' and 'Male'
#'
#' @export
splitGender <- function(.data) {
  assigned <- filter(.data, grepl("[^All]", gender))
  genders <- unique(assigned$gender)

  all <- filter(.data, grepl("All", gender))
  for(value in genders) {
    assigned <- rbind(assigned, mutate(all, gender = value))
  }

  assigned
#
#
#   female <- .data %>%
#     filter(grepl("(Female|All)", gender)) %>%
#     mutate(gender = "Female")
#
#   male <- .data %>%
#     filter(grepl("(Male|All)", gender)) %>%
#     mutate(gender = "Male")
#
#   rbind(female, male)
}

#' Filter calibrated forms from the rest
#' @export
filterCalibrated <- function(.data) {
  filter(.data, form == "Calibrated")
}

#' Filter well-defined forms from the rest
#' @export
filterFree <- function(.data) {
  filter(.data, form != "Calibrated")
}

#' Crush numbered variables
#'
#' All numbered variables (i.e. b1-b16) are crushed into a single variable.
#' Assumes numbered variables are al of numeric type.
#' New variable is named 'betas'.
#'
#' @export
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
#' @export
makeCalibratedFactories <- function(rr, pc, mm) {
  message("Building and calibrating constrained factories... ", appendLF = FALSE)

  .data <- inner_join(
    x = rr,
    y = inner_join(
      x = pc, y = mm, by = c("region", "year", "gender", "age_group")
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

  message("Done")

  .data
}



#' Factory for AAF computer factories: conditions with well-defined rel. risk
#' @export
makeFreeFactories <- function(.data) {
  message("Building unconstrained factories... ", appendLF = F)
  .data %<>%
    mutate(
      base_risk = pmap(
        list(im  = im, gender = gender, form = form, betas = betas),
        makeBaseRisk
      )
    ) %>%
    mutate(x1 = ifelse(grepl("5...[2ZR]", im), 50, 100)) %>%
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

  message("Done")

  .data
}

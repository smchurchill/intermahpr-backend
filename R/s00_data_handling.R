#### Verify and wrangle data ---------------------------------------------------

#' Format relative risk input data to our desired specifications
#'
#'@description
#'Formats all variable names to upper case and imputes imputable missing data.
#'Intended for a relative risk sheet as described in the intermahp user guide.
#'
#'expected variables:
#'  Gender (paired with PC, DH data)
#'  Function (FP, Spline, with valid IM, or Step, for HIV only)
#'  B1-B16 (only if Function == FP. need nonzero betas, rest can be missing)
#'  IM
#'  Condition
#'  Outcome
#'  RR_FD
#'  BingeF
#'
#'@param .data is a dataset that can be converted to tibble with the expected
#'  variables
#'
#'@export
#'

format_v1_rr <- function(.data) {
  .data <- tibble::as.tibble(.data)
  names(.data) <- do.call(stringr::str_to_upper, list(names(.data)))
  BETAS <- names(.data)[grep("[0-9]$", names(.data))]
  EXPECTED <- c(
    "IM", "CONDITION", "GENDER", "OUTCOME",
    "RR_FD", "BINGEF", "FUNCTION", BETAS
  )
  MISSING <- EXPECTED[!(EXPECTED %in% names(.data))]
  if(length(MISSING) > 0) {
    stop(
      "Missing variables from the supplied relative risk sheet: ",
      paste(MISSING, collapse = ", ")
    )
  }

  for(var in names(.data)) {
    .data[var][is.missing(.data[var])] <- impute_with(var)
  }

  readr::type_convert(.data[EXPECTED])
}

#' Format prevalence and consumption input data to our desired specifications
#'
#'@description
#'Formats all variable names to upper case and imputes imputable missing data.
#'Intended for a prevalence and consumption sheet as described in the intermahp
#'user guide.
#'
#' expected variables:
#'    P_FD (Proportion, Former Drinkers)
#'    P_BD (Proportion, Binge Drinkers)
#'    P_LA (Proportion, Lifetime Abstainers)
#'    P_CD (Proportion, Current Drinkers)
#'    Gamma Distribution Parameters:
#'      PCC_litres_year (Per Capita Consumption)
#'      Population
#'      Relative_Consumption
#'      Correction_Factor
#'    Gender (in order to pair with RR curves)
#'    Region
#'    Year
#'    Gender
#'    Age_Group
#'
#'@param .data is a dataset that can be converted to tibble with the expected
#'  variables
#'
#'@export
#'

format_v1_pc <- function(.data) {
  .data <- tibble::as.tibble(.data)
  names(.data) <- do.call(stringr::str_to_upper, list(names(.data)))

  EXPECTED <- c(
    "REGION", "YEAR", "GENDER", "AGE_GROUP", "POPULATION", "PCC_LITRES_YEAR",
    "CORRECTION_FACTOR", "RELATIVE_CONSUMPTION", "P_LA", "P_FD", "P_CD", "P_BD"
  )
  MISSING <- EXPECTED[!(EXPECTED %in% names(.data))]
  if(length(MISSING) > 0) {
    stop(
      "Missing variables from the supplied prevalence/consumption sheet: ",
      paste(MISSING, collapse = ", ")
    )
  }

  for(var in names(.data)) {
    .data[var][is.missing(.data[var])] <- impute_with(var)
  }

  readr::type_convert(.data[EXPECTED])
}

#' Format Count Data
#'
#'@description
#'Formats all variable names to upper case and imputes imputable missing data.
#'Intended for a count sheet
#'
#' expected variables:
#'    IM
#'    Region
#'    Year
#'    Gender
#'    Age_Group
#'    Outcome
#'    Count
#'
#'@param .data is a dataset that can be converted to tibble with the expected
#'  variables
#'
#'@importFrom magrittr %>% %<>%
#'
#'@export
#'

format_v1_dh <- function(.data) {
  .data <- tibble::as.tibble(.data)
  names(.data) <- do.call(stringr::str_to_upper, list(names(.data)))

  EXPECTED <- c(
    "IM", "CONDITION", "REGION", "YEAR", "GENDER", "AGE_GROUP",
    "OUTCOME", "COUNT"
  )
  MISSING <- EXPECTED[!(EXPECTED %in% names(.data))]
  if(length(MISSING) > 0) {
    stop(
      paste0(
        "Missing and needed variables from the supplied counts sheet: "
      ),
      paste(MISSING, collapse = ", ")
    )
  }

  for(var in names(.data)) {
    .data[var][is.missing(.data[var])] <- impute_with(var)
  }

  .data <- readr::type_convert(.data[EXPECTED])
  SUPPORTED <- paste0(
        "(7...1",
        "|3...[12]",
        "|1...[123]",
     "|[69]...[1-5]",
        "|8...[1-6]",
    "|[245]...[1-7])")
  IGNORE <- filter(.data, !grepl(SUPPORTED, IM))
  ignore_message <- paste0(
    "The following IMs are not supported by InterMAHP and will be ignored:\n",
    capture.output(unique(IGNORE$IM)),
    collapse = "\n"
  )
  warning(ignore_message)
  filter(.data, grepl(SUPPORTED, IM))
}

#### Derive params for initial AAF evaluation ----------------------------------

#' Derives relative risk curves from relative risk data
#'
#'@param .data tibble as formatted in format_v*_rr
#'@param ext logical indicator of extrapolation, (TRUE = linear, FALSE = capped)
#'
#'@importFrom magrittr %>% %<>%
#'@importFrom purrr map map2 pmap
#'
#'@export
#'

derive_v1_rr <- function(.data, ext) {
  NO_B <- .data[-grep("[0-9]", names(.data))]
  YES_B <- .data[grep("[0-9]$", names(.data))]

  BLIST <- split(as.matrix(YES_B), 1:nrow(YES_B))
  NO_B$BETAS <- BLIST

  NO_B %>%
    mutate(
      BASE_RR = pmap(
        list(IM, GENDER, FUNCTION, BETAS),
        base_rr_factory)
    ) %>%
    mutate(
      EXT = ifelse(
        FUNCTION == "Spline" & GENDER == "Female",
        FALSE, ext)
    ) %>%
    mutate(
      X1 = ifelse(
        grepl("5...2", IM),
        50, 100
      )
    ) %>%
    mutate(
      X2 = X1 + 50
    ) %>%
    mutate(
      Y1 = map2_dbl(BASE_RR, X1, ~.x(.y)),
      Y2 = map2_dbl(BASE_RR, X2, ~.x(.y))
    ) %>%
    mutate(
      CMP_SLOPE = (Y2-Y1)/(X2-X1)
    ) %>%
    mutate(
      SLOPE = ifelse(EXT, CMP_SLOPE, 0)
    ) %>%
    mutate(
      LNXT_RR = pmap(
        list(BASE_RR, X2, Y2, SLOPE),
        linear_extrapolation_factory
      )
    ) %>%
    mutate(
      BNGD_RR = pmap(
        list(IM, BINGEF, LNXT_RR),
        binge_risk_factory
      )
    )
}


#' Derives probability distribution parameters from Prevalence and Consumption
#' data
#'
#' From a given prevalence and consumption dataframe (or similar), computes re-
#' -quired data for AAF computation including Per Capita Consumption among
#' drinkers, and Gamma Shape/Scale parameters.  Also computes constants
#' necessary for AAF integration, namely the Deflation factor and Binge-Split
#' ratios.
#'
#' Gamma parameters are derivable entirely from prevalence and consumption,
#' no additional input is necessary.
#'
#'@param .data tibble as formatted in format_v*_pc
#'
#'@param bb Gender stratified binge barrier
#'@param lb Lower bound
#'@param ub Upper bound
#'@param gc Gender stratified gamma constants
#'
#' Constants used:
#'  0.002739726 = 1/365
#'  0.7893 = unit conversion (millilitres to grams) for ethanol
#'
#'@return A tibble with additional variables:
#'   PCC_G_DAY: Per capita consumption, converted from litres/year to grams/day
#'   DRINKERS: Total drinkers in region, computed as POPULATION*P_CD
#'   PCAD: Average consumption among drinkers
#'   PCC_AMONG_DRINKERS: Computed from PCAD and relative consumption
#'   GAMMA_CONSTANT: The distribution of drinkers at drinking level of x grams
#'     per day is fit by a gamma distribution whose mean and standard deviation
#'     are related by GAMMA_CONSTANT, stratified by gender
#'   GAMMA_SHAPE: Computed from GAMMA_CONSTANT
#'   GAMMA_SCALE: Computed from GAMMA_CONSTANT
#'   BB: Binge barrier
#'   LB: Lower bound of integration (0.03)
#'   UB: Upper bound of integration (user-specified, typical values are 150,250)
#'   NC: Normalizing constant for the gamma distribution (how much of the gamma
#'   distribution is between LB and UB?)
#'   DF: Deflation factor for the gamma distribution (DF = P_CD/NC)
#'   P_BAT: Proportion of bingers above threshold (BB)
#'   R1, R2: Ratios computed at this step and used in AAF computations
#'
#'@importFrom magrittr "%>%" "%<>%"
#'@importFrom dplyr group_by mutate
#'@importFrom tibble add_column
#'
#'@export
#'

derive_v1_pc <- function(.data, bb, lb, ub, gc) {
  .data %>%
    group_by(REGION, YEAR) %>%
    mutate(
      PCC_G_DAY = PCC_LITRES_YEAR * 1000 *
        0.7893 * CORRECTION_FACTOR * 0.002739726,
      DRINKERS = POPULATION * P_CD,
      PCAD = PCC_G_DAY * sum(POPULATION) / sum(DRINKERS),
      PCC_AMONG_DRINKERS = RELATIVE_CONSUMPTION * PCAD * sum(DRINKERS) /
        sum(RELATIVE_CONSUMPTION*DRINKERS)
      )%>%
    ungroup %>%
    mutate(
      GAMMA_CONSTANT = map_dbl(GENDER, ~`[[`(gc, .x))
    ) %>%
    mutate(
      GAMMA_SHAPE = 1/GAMMA_CONSTANT,
      GAMMA_SCALE = GAMMA_CONSTANT*PCC_AMONG_DRINKERS,
      BB = map_dbl(GENDER, ~`[[`(bb, .x)),
      LB = lb,
      UB = ub
    ) %>%
    mutate(
      GLB = pgamma(q = LB, shape = GAMMA_SHAPE, scale = GAMMA_SCALE),
      GBB = pgamma(q = BB, shape = GAMMA_SHAPE, scale = GAMMA_SCALE),
      GUB = pgamma(q = UB, shape = GAMMA_SHAPE, scale = GAMMA_SCALE)
    ) %>%
    mutate(
      NC = GUB - GLB
    ) %>%
    mutate(
      DF = P_CD / NC
    ) %>%
    mutate(
      N_GAMMA = pmap(
        list(GAMMA_SHAPE, GAMMA_SCALE, DF),
        normalized_gamma_factory
      ),
      P_BAT = DF * (GUB - GBB)
    ) %>%
    mutate(
      R1 = (P_CD - P_BD)  / (P_CD - P_BAT),
      R2 = (P_BD - P_BAT) / (P_CD - P_BAT)
    )
}


#' Derives count data for use in calibration and portioning
#'
#'@param dh Deaths/Hosps as returned by format_v1_dh
#'@param pc Prev/Cons as returned by derive_v1_pc
#'
#'@importFrom magrittr %>% %<>%
#'
#'@export

derive_v1_dh <- function(dh, pc) {
  PC <- pc[
    c("REGION", "YEAR", "GENDER", "AGE_GROUP", "DRINKERS",
      "BB", "LB", "UB", "N_GAMMA")]
  DH <- dh %>%
  dplyr::inner_join(
    PC,
    by = c("REGION", "YEAR", "GENDER", "AGE_GROUP")
  )
  DH
}


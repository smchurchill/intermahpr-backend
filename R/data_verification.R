impute_with <- function(var) {
  if(grepl("B[0-9]", var)) {var = "BETA"}
  switch(var,
         IM = "(0).(0)",
         CONDITION = "Unspecified",
         REGION = "Unspecified",
         YEAR = "Unspecified",
         AGE_GROUP = "Unspecified",
         OUTCOME = "Combined",
         FUNCTION = "FP",
         RR_FD = 1,
         BINGEF = 1,
         BETA = 0)
}

is.missing <- function(obs) {
  is.na(obs) | is.null(obs) | obs == "."
}

impute_missing <- function(var_obs) {
  TF <- is.missing(var_obs)
  var_obs[TF] <- impute_with(var_obs)
}

format_rr <- function(RR_in) {
  RR <- tibble::as.tibble(RR_in)
  names(RR) <- do.call(stringr::str_to_upper, list(names(RR)))
  expected_var <- c("IM", "CONDITION", "GENDER", "OUTCOME", "RR_FD", "BINGEF",
                    "FUNCTION", do.call(paste0, list(rep("B", 16), 1:16)))
  missing_var <- expected_var[!(expected_var %in% names(RR))]

  RR[, missing_var] <- NA

  for(var in names(RR)) {
    RR[var][is.missing(RR[var])] <- impute_with(var)
  }

  readr::type_convert(RR)
}

format_v0_pc <- function(PC_in) {
  PC <- tibble::as.tibble(PC_in)
  names(PC) <- do.call(stringr::str_to_upper, list(names(PC)))
  ## expected variables:
  ##  Needed:
  ##    Proportion, Former Drinkers
  ##    Proportion, Binge Drinkers
  ##    Proportion, Lifetime Abstainers
  ##    Proportion, Current Drinkers
  ##    Gamma Distribution Parameters:
  ##      Per Capita Consumption
  ##      Population
  ##      Relative Consumption
  ##      Correction Factor
  ##    Gender (in order to pair with RR curves)
  ##  Imputable as "Unspecified":
  ##    Region
  ##    Year
  ##    Gender
  ##    Age Group
  ##

  expected_var <- c("REGION", "YEAR", "GENDER", "AGE_GROUP", "POPULATION",
                    "PCC_LITRES_YEAR", "CORRECTION_FACTOR",
                    "RELATIVE_CONSUMPTION", "P_LA", "P_FD", "P_CD", "P_BD")
  missing_var <- expected_var[!(expected_var %in% names(PC))]
  needed_var <- expected_var[c(3,5:12)]
  missing_and_needed <- intersect(missing_var, needed_var)
  if(length(missing_and_needed) > 0) {
    stop(paste0("Missing and needed variables from the supplied prevalence/",
                "consumption sheet: "),
                paste(missing_and_needed, collapse = ", "))
  }

  PC[, missing_var] <- NA

  for(var in names(PC)) {
    PC[var][is.missing(PC[var])] <- impute_with(var)
  }

  readr::type_convert(PC)
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
#'@param pc_input is assumed type data.table with variable list of c(Year,
#' Region, Gender, Age_group, Population, PCC_litres_year, Correctio\
#' n_factor, Relative_consumption, P_LA, P_FD, P_CD, P_BD).
#' Variables added to this list in the returned data.table are c(PCC_among_dri\
#' nkers, Gamma_shape, Gamma_scale, nc, df, p_bat, R1, R2).  This assumption is
#' safe because this function should be internally called only, not exposed.
#'
#' Note that the bundled data set pc_default satisfies these constraints.
#'
#' 0.002739726 = 1/365
#' 0.7893 = unit conversion (millilitres to grams) for ethanol
#'
#'@return A data.table with columns Year, Region, Gender, Age_group, Population,
#' PCC_litres_year, Correction_factor, Relative_consumption, P_LA, P_FD, P_CD,
#' P_BD, PCC_among_drinkers, Gamma_shape, Gamma_scale, nc, df, p_bat, R1, R2
#'
#'
#'


derive_v0_pc <- function(PC_in,
                         bb = c(53.8, 67.25),
                         lb = c(0.03, 0.03),
                         ub = c(250,  250)) {
  PC <- PC_in %>%
    group_by(REGION, YEAR) %>%
    mutate(PCC_G_DAY = PCC_LITRES_YEAR * 1000 *
             0.7893 * CORRECTION_FACTOR * 0.002739726,
           DRINKERS = POPULATION * P_CD,
           PCAD = PCC_G_DAY * sum(POPULATION) / sum(DRINKERS),
           PCC_AMONG_DRINKERS = RELATIVE_CONSUMPTION * PCAD * sum(DRINKERS) /
             sum(RELATIVE_CONSUMPTION*DRINKERS),
           GAMMA_CONSTANT = if_else(GENDER == "FEMALE", 1.258^2, 1.171^2),
           GAMMA_SHAPE = 1/GAMMA_CONSTANT,
           GAMMA_SCALE = GAMMA_CONSTANT*PCC_AMONG_DRINKERS,
           INDEX = if_else(GENDER == "FEMALE", 1, 2),
           BB = bb[INDEX],
           LB = lb[INDEX],
           UB = ub[INDEX],
           NC = integrate(function(x) dgamma(x,
                                             shape = GAMMA_SHAPE,
                                             scale = GAMMA_SCALE),
                          lower = LB,
                          upper = UB)$value,
           DF = P_CD / NC,
           P_BAT = integrate(function(x) DF * dgamma(x,
                                                     shape = GAMMA_SHAPE,
                                                     scale = GAMMA_SCALE),
                             lower = BB,
                             upper = UB)$value,
           R1 = (P_CD - P_BD)  / (P_CD - P_BAT),
           R2 = (P_BD - P_BAT) / (P_CD - P_BAT))
  PC
}

#' Default Prevalence and Consumption Data
#'
#' Default prevalence and consumption data from BC and Canada, 2015.input to
#' generate relative risk curves.  Standard formatting is described in the
#' InterMAHP user guides
#'
#' @docType data
#'
#' @usage data(pcr_default)
#'
#'
"pc_default"


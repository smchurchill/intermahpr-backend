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
#'@import data.table
#'
#'

derive_params_from_pc <- function(pc = intermahpr::pc_default,
                                  bb = c(53.8, 67.25),
                                  lb = c(0.03, 0.03),
                                  ub = c(250,  250)){
  DT <- data.table::data.table(pc)

# Obtain mean consumption among drinkers in grams/day
  # Convert PCC from litres/year to g/day, adjust based on correction factor
  DT[, PCC_g_day := PCC_litres_year*1000*0.7893*Correction_factor*0.002739726]

  # Number of drinkers in each age group
  DT[, drinkers := Population * P_CD]

  # Alcohol consumption by age group
  DT[, pcad :=
         PCC_g_day * sum(Population) /
         sum(drinkers), by = c("Region", "Year")]

  # Mean consumption by age group
  DT[, PCC_among_drinkers :=
         Relative_consumption*pcad*sum(drinkers) /
         sum(Relative_consumption*drinkers), by = c("Region", "Year")]

# Use mean consumption to produce gamma shape and scale parameters
  # First, the by-gender linear relationship between sigma**2 and mu**2
  fgc = 1.258^2 # female gamma constant
  rfgc = 1/fgc  # recioprocal ...
  mgc = 1.171^2 # male ...
  rmgc = 1/mgc  # ...
  # Used to compute Gamma shape and scale parameters
  DT[Gender == "Female",      `:=`(Gamma_shape = rfgc,
                                   Gamma_scale = fgc*PCC_among_drinkers,
                                   BB = bb[1],
                                   LB = lb[1],
                                   UB = ub[1])
   ][Gender == "Male",        `:=`(Gamma_shape = rmgc,
                                   Gamma_scale = mgc*PCC_among_drinkers,
                                   BB = bb[2],
                                   LB = lb[2],
                                   UB = ub[2])]
  # Then compute deflation factors and binge-splits for final integrands
  DT[, `:=` (nc = mapply(function(y,z) {
                   integrate(function(x) dgamma(x, shape = y, scale = z),
                     lower = intermahpr_constants$lower_bound,
                     upper = intermahpr_constants$upper_bound)$value
                   },
                   Gamma_shape,
                   Gamma_scale)
            )
    ]

  DT[, df := P_CD / nc]

  DT[, `:=`(p_bat = mapply(function(shape, scale, deflate, binge) {
                      integrate(function(x) dgamma(x,
                                                   shape = shape,
                                                   scale = scale),
                        lower = binge,
                        upper = intermahpr_constants$upper_bound)$value
                    },
                    Gamma_shape,
                    Gamma_scale,
                    df,
                    bb)
           )
    ]

  DT[, `:=`(R1 = (P_CD - P_BD)  / (P_CD - p_bat),
            R2 = (P_BD - p_bat) / (P_CD - p_bat))]

# Remove extraneous columns
  DT[, pcad := NULL][, drinkers := NULL]

  DT
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


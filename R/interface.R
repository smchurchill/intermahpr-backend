#' Collect and assemble AAF data from formatted RR and PC data
#'
#'@param pc  Prevalence / Consumption input
#'@param rr  Relative Risk input
#'@param ext logical, extrapolate linearly?
#'@param lb  Double, consumption lower bound
#'@param ub  Double, consumption upper bound
#'@param bb  Double vector, Binge consumption level, Gender stratified
#'@param gc  Gamma constant.  The linear relationship between mean and standard
#'  deviation within the gamma distribution that describes consumption among
#'  current drinkers.  Stratified by gender -- names(gc) must match levels of
#'  rr$GENDER and pc$GENDER.
#'
#'

assemble <- function(pc, rr, ext, lb, ub, bb, gc) {
  RRD <- derive_v0_rr(rr = rr, ext = ext)
  PCD <- derive_v0_pc(pc = pc, bb = bb, lb = lb, ub = ub, gc = gc)

  join_pc_rr(pc = PCD, rr = RRD)
}

#' Compute AAFs for all conditions as separated by the given cutpoints
#'
#'@param aaf_table A tibble as returned by assemble
#'@param cuts A list of double vectors indexed by aaf_table$GENDER
#'
#'@importFrom purrr pmap map map2
#'@importFrom magrittr "%>%" "%<>%"
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

#' Takes all possible inputs, computes specified AAFs, returns TMI
#'
#'@description
#'The most granular AAF computation.  All inputs have default values, and levels
#'of Gender are inspected before any computations are done.
#'
#'@param pc  Prevalence / Consumption input
#'@param rr  Relative Risk input
#'@param ext logical.  Answer to the question: "extrapolate linearly?"
#'@param lb  Double, consumption lower bound, given in
#'  units of grams ethanol per day
#'@param ub  Double, consumption upper bound, given in
#'  units of grams ethanol per day
#'@param bb  Double vector, Binge consumption level, Gender stratified, given in
#'  units of grams ethanol per day. Stratified by gender -- names(bb) must match
#'  levels of rr$GENDER and pc$GENDER.
#'@param gc  Gamma constant. The linear relationship between mean and standard
#'  deviation within the gamma distribution that describes consumption among
#'  current drinkers. Stratified by gender -- names(gc) must match levels of
#'  rr$GENDER and pc$GENDER.
#'@param cb  Consumption barriers. These define light/moderate/heavy drinking
#'  groups. Given in units of grams ethanol per day. Stratified by gender --
#'  names(cb) must match levels of rr$GENDER and pc$GENDER.
#'
#'@export
#'

intermahpr_raw <- function(pc = intermahpr::pc_default,
                           rr = intermahpr::rr_default,
                           ext = TRUE, lb = 0.03, ub = 250,
                           bb = list("Female" = 53.8, "Male" = 67.25),
                           gc = list("Female" = 1.258^2, "Male" = 1.171^2),
                           cb = list("Female" = c(13.45, 26.9),
                                     "Male" = c(20.2, 40.4))) {
  PCF <- format_v0_pc(pc = pc)
  RRF <- format_v0_rr(rr = rr)
  PC_LEVELS <- levels(as.factor(PCF$GENDER))
  RR_LEVELS <- levels(as.factor(RRF$GENDER))
  BB_LEVELS <- levels(as.factor(names(bb)))
  GC_LEVELS <- levels(as.factor(names(gc)))
  CB_LEVELS <- levels(as.factor(names(cb)))
  LEVEL_MATCHES <- sapply(list(RR_LEVELS, BB_LEVELS, GC_LEVELS, CB_LEVELS),
                          FUN = identical, PC_LEVELS)
  LEVELS_IDENTICAL <- all(LEVEL_MATCHES)

  if(!LEVELS_IDENTICAL) {
    LEVELS_MESSAGE <- paste0("Prevalence/Consumption Gender levels: ",
                             capture.output(PC_LEVELS),
                             "\n",
                             "Relative Risk Gender levels:          ",
                             capture.output(RR_LEVELS),
                             "\n",
                             "Binge Barrier Gender levels:          ",
                             capture.output(BB_LEVELS),
                             "\n",
                             "Gamma Constant Gender levels:         ",
                             capture.output(GC_LEVELS),
                             "\n",
                             "Consumption Barrier Gender levels:    ",
                             capture.output(CB_LEVELS),
                             collapse = "\n")

    stop(paste0("Levels of Gender variable within Prevalence/Consumption, R",
                   "elative Risk, Binge Barrier, Gamma Constant, and Consumpti",
                   "on Barriers must be equal.  The levels input are:\n",
                   collapse = "\n"),
            LEVELS_MESSAGE)
  }

  ASSEMBLED <- assemble(pc = PCF, rr = RRF,
                        ext = ext, lb = lb, ub = ub,
                        bb = bb, gc = gc)
  COMPUTED <- compute_aafs(aaf_table = ASSEMBLED, cuts = cb)
  COMPUTED
}

#' Helper function that tidies base AAF into a certain outcome type
#'
#'
#'@importFrom magrittr "%>%" "%<>%"
#'@importFrom dplyr filter
#'


outcome_splitter <- function(aaf_table, outcome) {
  aaf_table %<>%
    filter(OUTCOME == outcome | OUTCOME == "Combined") %>%
    mutate(AAF_LD = sapply(AAF_GRP, `[[`, 1),
           AAF_MD = sapply(AAF_GRP, `[[`, 2),
           AAF_HD = sapply(AAF_GRP, `[[`, 3)) %>%
    select(REGION, YEAR, GENDER, AGE_GROUP, IM, CONDITION,
           AAF_FD, AAF_LD, AAF_MD, AAF_HD, AAF_TOTAL)

  aaf_table
}

#' Integrate that takes as input a list funs of functions, a vector vlower of
#' "lower" values and a vector vupper of "upper" values, assumes that
#' length(vlower) = length(vupper) = length(funs), and returns a
#' vector of integrals of funs between the lower and upper values.
#'
#'
#'

vintegrate <- function(funs, vlower, vupper) {
  llower <- length(vlower)
  lupper <- length(vupper)
  lfuns <- length(funs)

  values <- rep(-1, lfuns)

  for(i in 1:lfuns) {
    values[[i]] <- integrate(funs[[i]], vlower[[i]], vupper[[i]])$value
  }

  values
}

#' Helper function that tidies base AAF into the relevant prevalence and
#' consumption data
#'
#'@importFrom dplyr distinct select case_when
#'

extract_prevcons <- function(aaf_table) {
  aaf_table %<>%
    mutate(LM = sapply(CUTS, `[[`, 2),
           MH = sapply(CUTS, `[[`, 3)) %>%
    distinct(REGION, YEAR, GENDER, AGE_GROUP, POPULATION, PCC_AMONG_DRINKERS,
             GAMMA_SHAPE, GAMMA_SCALE, P_LA, P_FD, P_CD, EXT, LB, LM, MH,
             UB, .keep_all = TRUE) %>%
    select(REGION, YEAR, GENDER, AGE_GROUP, POPULATION, PCC_AMONG_DRINKERS,
           GAMMA_SHAPE, GAMMA_SCALE, P_LA, P_FD, P_CD, EXT, LB, LM, MH, UB,
           N_GAMMA)

  aaf_table %<>%
    mutate(EXTRAPOLATION = case_when(EXT == TRUE  ~ "linear",
                                     EXT == FALSE ~ "capped")) %>%
    add_column(P_LD = vintegrate(aaf_table[["N_GAMMA"]],
                                 aaf_table[["LB"]],
                                 aaf_table[["LM"]]),
               P_MD = vintegrate(aaf_table[["N_GAMMA"]],
                                 aaf_table[["LM"]],
                                 aaf_table[["MH"]]),
               P_HD = vintegrate(aaf_table[["N_GAMMA"]],
                                 aaf_table[["MH"]],
                                 aaf_table[["UB"]]))%>%
    mutate(P_CD_SUM = P_LD + P_MD + P_HD) %>%
    select(REGION, YEAR, GENDER, AGE_GROUP, POPULATION, PCC_AMONG_DRINKERS,
           GAMMA_SHAPE, GAMMA_SCALE, P_LA, P_FD, P_CD, P_LD, P_MD, P_HD,
           EXTRAPOLATION, P_CD_SUM)

  aaf_table
}


#' Takes Base InterMAHP inputs, returns formatted data
#'
#'@param RelativeRisks rr
#'@param PrevalenceConsumption pc
#'@param FemaleLightModerateBarrier flm
#'@param FemaleModerateHeavyBarrier fmh
#'@param FemaleBingeBarrier fbb
#'@param MaleLightModerateBarrier mlm
#'@param MaleModerateHeavyBarrier mmh
#'@param MaleBingeBarrier mbb
#'@param UpperBound ub
#'@param Extrapolation ext
#'@param OutputPath op
#'
#'@export
#'

intermahpr_base <- function(RelativeRisks = intermahpr::rr_default,
                            PrevalenceConsumption,
                            FemaleLightModerateBarrier,
                            FemaleModerateHeavyBarrier,
                            FemaleBingeBarrier,
                            MaleLightModerateBarrier,
                            MaleModerateHeavyBarrier,
                            MaleBingeBarrier,
                            UpperBound,
                            Extrapolation,
                            OutputPath = NULL,
                            FilePrefix = "") {

  BB <- list("Female" = FemaleBingeBarrier, "Male" = MaleBingeBarrier)
  GC <- list("Female" = 1.258^2, "Male" = 1.171^2)
  CB <- list("Female" = c(FemaleLightModerateBarrier,
                          FemaleModerateHeavyBarrier),
             "Male" = c(MaleLightModerateBarrier,
                        MaleModerateHeavyBarrier))
  AAF_OUT <- intermahpr_raw(rr = RelativeRisks, pc = PrevalenceConsumption,
                            ext = Extrapolation, lb = 0.03, ub = UpperBound,
                            bb = BB, gc = GC, cb = CB)
  InterMAHP_AAFs_morbidity <- outcome_splitter(AAF_OUT, "Morbidity")
  InterMAHP_AAFs_mortality <- outcome_splitter(AAF_OUT, "Mortality")
  InterMAHP_prev_cons_output <- extract_prevcons(AAF_OUT)

  if(file.exists(OutputPath)) {
    readr::write_csv(x = InterMAHP_AAFs_morbidity,
                     path = file.path(OutputPath,
                                      paste0(FilePrefix,
                                             "InterMAHP_AAFs_morbidity.csv")))
    readr::write_csv(x = InterMAHP_AAFs_mortality,
                     path = file.path(OutputPath,
                                      paste0(FilePrefix,
                                             "InterMAHP_AAFs_mortality.csv")))
    readr::write_csv(x = InterMAHP_AAFs_morbidity,
                     path = file.path(OutputPath,
                                      paste0(FilePrefix,
                                             "InterMAHP_prev_cons_output.csv")))
  }
  else if(!is.null(OutputPath)) {
    message(
      paste0(
        "Path does not exist.  Returning list of output tables instead. If you",
        " haven't explicitly assigned output to a variable, this list is store",
        "d in the variable '.Last.value'"
      )
    )
  }

  invisible(
    list(
      "InterMAHP_AAFs_morbidity" = InterMAHP_AAFs_morbidity,
      "InterMAHP_AAFs_mortality" = InterMAHP_AAFs_mortality,
      "InterMAHP_prev_cons_output" = InterMAHP_prev_cons_output
    )
  )
}

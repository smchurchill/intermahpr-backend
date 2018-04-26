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

intermahpr_raw <- function(
  pc = intermahpr::pc_default,
  rr = intermahpr::rr_default,
  ext = TRUE, lb = 0.03, ub = 250,
  bb = list("Female" = 53.8, "Male" = 67.25),
  gc = list("Female" = 1.258^2, "Male" = 1.171^2),
  cb = list("Female" = c(13.45, 26.9), "Male" = c(20.2, 40.4))
){
  PCF <- format_v0_pc(pc = pc)
  RRF <- format_v0_rr(rr = rr)
  PC_LEVELS <- levels(as.factor(PCF$GENDER))
  RR_LEVELS <- levels(as.factor(RRF$GENDER))
  BB_LEVELS <- levels(as.factor(names(bb)))
  GC_LEVELS <- levels(as.factor(names(gc)))
  CB_LEVELS <- levels(as.factor(names(cb)))
  LEVEL_MATCHES <- sapply(
    X = list(RR_LEVELS, BB_LEVELS, GC_LEVELS, CB_LEVELS),
    FUN = identical, PC_LEVELS
  )
  LEVELS_IDENTICAL <- all(LEVEL_MATCHES)

  if(!LEVELS_IDENTICAL) {
    LEVELS_MESSAGE <- paste0(
      "Prevalence/Consumption Gender levels: ",
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
      collapse = "\n"
    )
    stop(
      paste0(
        "Levels of Gender variable within Prevalence/Consumption, Relative ",
        "Risk, Binge Barrier, Gamma Constant, and Consumption Barriers must ",
        "be equal.  The levels input are:\n",
        collapse = "\n"),
      LEVELS_MESSAGE
    )
  }

  ASSEMBLED <- assemble(
    pc = PCF, rr = RRF, ext = ext, lb = lb, ub = ub, bb = bb, gc = gc
  )
  COMPUTED <- compute_aafs(aaf_table = ASSEMBLED, cuts = cb)
  COMPUTED
}

#' Takes Base InterMAHP inputs, returns formatted data
#'
#'@description
#'  Basic interface with intermahpr package.  Documentation TBW.
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

intermahpr_base <- function(
  RelativeRisks = intermahpr::rr_default,
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
  FilePrefix = ""
){
  BB <- list("Female" = FemaleBingeBarrier, "Male" = MaleBingeBarrier)
  GC <- list("Female" = 1.258^2, "Male" = 1.171^2)
  CB <- list(
    "Female" = c(
      FemaleLightModerateBarrier,
      FemaleModerateHeavyBarrier
    ),
    "Male" = c(
      MaleLightModerateBarrier,
      MaleModerateHeavyBarrier
    )
  )
  AAF_OUT <- intermahpr_raw(
    rr = RelativeRisks, pc = PrevalenceConsumption, ext = Extrapolation,
    lb = 0.03, ub = UpperBound, bb = BB, gc = GC, cb = CB
  )

  SPLIT <- outcome_splitter(AAF_OUT)
  NAMED <- name_cuts(SPLIT)

  InterMAHP_AAFs_morbidity <- NAMED[NAMED$OUTCOME == "Morbidity", ]
  InterMAHP_AAFs_mortality <- NAMED[NAMED$OUTCOME == "Mortality", ]
  InterMAHP_prev_cons_output <- extract_prevcons(AAF_OUT)

  output <- list(
    "InterMAHP_AAFs_morbidity" = InterMAHP_AAFs_morbidity,
    "InterMAHP_AAFs_mortality" = InterMAHP_AAFs_mortality,
    "InterMAHP_prev_cons_output" = InterMAHP_prev_cons_output
  )


  if(is.null(OutputPath)) {return(output)}
  else if(file.access(file.path(OutputPath), 2) == 0) {
    readr::write_csv(
      x = InterMAHP_AAFs_morbidity,
      path = file.path(
        OutputPath,
        paste0(
          FilePrefix,
          "InterMAHP_AAFs_morbidity.csv"
        )
      )
    )
    readr::write_csv(
      x = InterMAHP_AAFs_mortality,
      path = file.path(
        OutputPath,
        paste0(
          FilePrefix,
          "InterMAHP_AAFs_mortality.csv"
        )
      )
    )
    readr::write_csv(
      x = InterMAHP_AAFs_morbidity,
      path = file.path(
        OutputPath,
        paste0(
          FilePrefix,
          "InterMAHP_prev_cons_output.csv"
        )
      )
    )
  }
  else {
    message(
      paste0(
        "Path does not exist, or R does not have write permission to the speci",
        "fied path. Returning list of output tables instead. If you",
        " haven't explicitly assigned output to a variable, this list is store",
        "d in the variable '.Last.value'"
      )
    )
  }

  output
}

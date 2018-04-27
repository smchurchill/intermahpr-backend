library(tidyverse)
library(magrittr)

dh_in <- readr::read_csv(file.path("data-raw", "dh.csv"),
                         col_types = "??????dddd????????")

dh <- dh_in %>%
  rename(region = Province) %>%
  rename(condition = Condition_Alcohol)

dh$Outcome <- "Morbidity"

dh

impc <- readr::read_csv(file.path("data-raw", "impc.csv"))

PC <- impc %>%
  intermahpr::format_v0_pc() %>%
  intermahpr::derive_v0_pc(
    bb = list("Female" = 50, "Male" = 60),
    lb = 0.03,
    ub = 250,
    gc = list("Female" = 1.258^2, "Male" = 1.171^2))

PC[c("REGION", "YEAR", "GENDER", "AGE_GROUP", "DRINKERS", "BB", "LB", "UB", "N_GAMMA")]

DH <- dh %>%
  intermahpr::format_v0_dh()

DH

DH %<>%
  intermahpr::derive_v0_dh(PC)

View(DH)


aaf <- intermahpr_base(RelativeRisks = rr_default, PrevalenceConsumption = pc_default, FemaleLightModerateBarrier = 10, FemaleModerateHeavyBarrier = 20, FemaleBingeBarrier = 53, MaleLightModerateBarrier = 15, MaleModerateHeavyBarrier = 30, MaleBingeBarrier = 65, UpperBound = 250, Extrapolation = TRUE)


aaf

aaf <- intermahpr_raw() %>%
  outcome_splitter() %>%
  name_cuts()

aaf

epi_sim <- epilepsy <- filter(aaf, IM == "(4).(5)")
epi_aaf <- epilepsy[c("AAF_LD", "AAF_MD", "AAF_HD", "AAF_TOTAL")]

epi_xaaf <- epi_aaf / epi_aaf$AAF_TOTAL
epi_sim[c("AAF_LD", "AAF_MD", "AAF_HD", "AAF_TOTAL")] <- epi_xaaf

epilepsy
epi_sim

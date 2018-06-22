library(tidyverse)
library(magrittr)

rr <- readr::read_csv("data-raw/rr_master.csv")
rr_p <- prepareRR(rr, T)

pc <- readr::read_csv("data-raw/pc_master.csv")
pc_p <- preparePC(pc, bb = list("Female" = 10, "Male" = 15))

dh <- readr::read_csv("data-raw/dh_master.csv")
dh_p <- prepareDH(dh)

shinymodel <- makeNewModel(rr = rr_p, pc = pc_p, dh = dh_p)

af <- formatForShinyOutput(shinymodel$scenarios$base)

names(shinymodel$scenarios)

scenario_analysis <- makeNewModel(rr = rr_p, pc = pc_p, dh = dh_p) %>%
  makeScenarios(scenario_names = c("MUP1", "MUP2"), scales = c(0.97, 0.95)) %>%
  distillModel()




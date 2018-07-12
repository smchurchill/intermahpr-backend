library(tidyverse)
library(magrittr)

rr <- readr::read_csv("C:/Users/samuelch.UVIC/Documents/shiny-inputs/rr_master.csv")
rr_p <- prepareRR(rr, T)

pc <- readr::read_csv("C:/Users/samuelch.UVIC/Documents/shiny-inputs/all_pc.csv")
pc_p <- preparePC(pc, bb = list("Female" = 53.8, "Male" = 67.25))
myStats <- dplyr::filter(pc_p, region == "20:CA" & year == 2014 & gender == "Male" & age_group == "35-64")
View(myStats)

myStats$gamma_shape
myStats$gamma_scale

pc <- readr::read_csv("C:/Users/samuelch.UVIC/Documents/shiny-inputs/pc_master.csv")
pc_p <- preparePC(pc, bb = list("Female" = 10, "Male" = 15))

dh <- readr::read_csv("C:/Users/samuelch.UVIC/Documents/shiny-inputs/all_dh.csv")

dh <- readr::read_csv("C:/Users/samuelch.UVIC/Documents/shiny-inputs/nl_dh.csv")

dh <- readr::read_csv("C:/Users/samuelch.UVIC/Documents/shiny-inputs/dh_master.csv")
dh_p <- prepareDH(dh)

shinymodel <- makeNewModel(rr = rr_p, pc = pc_p, dh = dh_p)

test2 <- shinymodel$scenarios$base %>% formatForShinyOutput()

test2 %>% filter(grepl("5...3", im))

test <- filter(shinymodel$scenarios$base, grepl("5...3", im))

test %<>% addCurrentFraction()

test

af <- formatForShinyOutput(shinymodel$scenarios$base)

names(shinymodel$scenarios)

scenario_analysis <- makeNewModel(rr = rr_p, pc = pc_p, dh = dh_p) %>%
  makeScenarios(scenario_names = c("MUP1", "MUP2"), scales = c(0.97, 0.95)) %>%
  distillModel()




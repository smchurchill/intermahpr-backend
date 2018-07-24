library(tidyverse)
library(magrittr)

rr <- readr::read_csv("C:/Users/samuelch.UVIC/Documents/shiny-inputs/rr_master.csv")
rr_p <- prepareRR(rr, T)

pc <- readr::read_csv("C:/Users/samuelch.UVIC/Documents/shiny-inputs/pc_master.csv")
pc_p <- preparePC(pc, bb = list("Female" = 10, "Male" = 15))

mm <- readr::read_csv("C:/Users/samuelch.UVIC/Documents/shiny-inputs/mm_master.csv")
mm_p <- preparemm(mm)

shinymodel <- makeNewModel(rr = rr_p, pc = pc_p, mm = mm_p)

shinymodel <- makeScenario(shinymodel, "Base", 1)

base_scen <- shinymodel$scenarios$Base

b2frac <- addGenderStratifiedIntervalFraction(
  base_scen,
  lower_strata = list("Female" = 10, "Male" = 20),
  upper_strata = list("Female" = 25, "Male" = 40)
)

test2 <- shinymodel$scenarios$base %>% formatForShinyOutput()

test2 %>% filter(grepl("5...3", im))

test <- filter(shinymodel$scenarios$base, grepl("5...3", im))

test %<>% addCurrentFraction()

test

af <- formatForShinyOutput(shinymodel$scenarios$base)

names(shinymodel$scenarios)

scenario_analysis <- makeNewModel(rr = rr_p, pc = pc_p, mm = mm_p) %>%
  makeScenarios(scenario_names = c("MUP1", "MUP2"), scales = c(0.97, 0.95)) %>%
  distillModel()


tryCatch(
  {
    print("first")
    stop("second")
    print("third")
    stop("fourth")
  }, error = function(e) {e$message}
)


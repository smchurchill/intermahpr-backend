library(tidyverse)
library(magrittr)

rr <- readr::read_csv("C:/Users/samuelch.UVIC/Documents/shiny-inputs/rr_master.csv")
pc <- readr::read_csv("C:/Users/samuelch.UVIC/Documents/shiny-inputs/pc_master.csv")
mm <- readr::read_csv("C:/Users/samuelch.UVIC/Documents/shiny-inputs/mm_master.csv")

pc_nogc <- pc
pc_nogc$gamma_constant <- NULL

rr_p <- prepareRR(rr, T)
pc_p <- preparePC(pc, bb = list("Female" = 10, "Male" = 15))
mm_p <- prepareMM(mm)

shinymodel <- makeNewModel(rr = rr_p, pc = pc_p, mm = mm_p) %>%
  makeScenario("Base", 1) %>%
  makeScenario("-2% consumption", .98)


ddd <- distillModel(shinymodel)


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


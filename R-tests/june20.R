library(tidyverse)
library(magrittr)

rr <- rr_master
rr_p <- prepareRR(rr, T)

pc <- pc_default
pc_p <- preparePC(pc)

dh <- readr::read_csv("data-raw/dh.csv", col_types = "??????dddd????????") %>%
  mutate(outcome = "Morbidity") %>%
  mutate(region = substring(Province, 4))
dh_p <- prepareDH(dh)

scenario_analysis <- makeNewModel(rr = rr_p, pc = pc_p, dh = dh_p) %>%
  makeScenarios(scenario_names = c("MUP1", "MUP2"), scales = c(0.97, 0.95)) %>%
  distillModel()

model <-

model$scenarios

MUPs <- makeScenarios(
  model,
  scenario_names = c("MUP1", "MUP2"),
  scales = c(0.97, 0.95)
)

dist <- distillModel(MUPs)

s1 <- MUPs$scenarios$MUP1

baseSc <- MUPs$scenarios$base

baseSc <- addTotalFraction(baseSc)

baseSc

basefd <- mutate(baseSc, aaf_fd = computeFormerFraction(baseSc))

computeFormerFraction(baseSc)
computeIntervalFraction(baseSc)
computeCurrentFraction(baseSc)
computeTotalFraction(baseSc)

computeCurrentFraction(s1)

names(model$scenarios)

map_dbl(MUPs$scenarios$base$current_fraction, ~.x(Inf)) %>% round(digits = 4)

for(i in 1:960) {
  print(i)
  print(object.size(model$model$current_fraction_factory[[i]], units = "Mb"))
  print(object.size(model$scenarios$base$current_fraction[[i]], units = "Mb"))
}


new <- generateScenario(model = scenarios, scenario_name = "what", scale = 1)

another <- generateScenario(model = new, scale = 0.97)

another$scenarios$what$current_fraction_factory


new




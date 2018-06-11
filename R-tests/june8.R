library(tidyverse)
library(magrittr)

rr <- rr_master
rr_p <- prepareRR(rr, T)

pc <- pc_default
pc_p <- preparePC(pc)

fn <- function() {
  print(pc)
}

fn()

dh <- readr::read_csv("data-raw/dh.csv", col_types = "??????dddd????????") %>%
  mutate(outcome = "Morbidity") %>%
  mutate(region = substring(Province, 4))
dh_p <- prepareDH(dh)

model <- makeNewModel(rr = rr_p, pc = pc_p, dh = dh_p)

model$scenarios$base

MUPs <- generateScenarios(
  model,
  scenario_names = c("MUP1", "MUP2"),
  scales = c(0.97, 0.95)
)

for(i in 1:960) {
  print(i)
  print(object.size(model$model$current_fraction_factory[[i]]) / 1000)
  print(object.size(model$scenarios$base$current_fraction[[i]]) / 1000)
}


new <- generateScenario(model = scenarios, scenario_name = "what", scale = 1)

another <- generateScenario(model = new, scale = 0.97)

another$scenarios$what$current_fraction_factory


new


scen1 <- scen0 %>%
  left_join(pc_p, by = c("region", "year", "gender", "age_group")) %>%
  mutate(
    fargs = pmap(
      list(p_fd = p_fd, mass = n_gamma, non_bingers = non_bingers, bingers = bingers, lb = lb, bb = bb, ub = ub),
      list
    )
  ) %>%
  mutate(
    current_fraction = map2(
      fargs,
      current_fraction_factory,
      ~.y(.x)
    )
  ) %>%
  mutate(
    former_fraction = map2(
      fargs,
      former_fraction_factory,
      ~.y(.x)
    )
  ) %>%
  mutate(aaf_cd = map2_dbl(ub, current_fraction, ~.y(.x)))

for(i in nrow(scen1):1) {
  print(i)
  print(scen1[i,]$current_fraction[[1]](250))
}


r_cal <- filterCalibrated(rr_p)
r_free <- filterFree(rr_p)
rr_free_factories <- makeFreeFactories(r_free)




args <- pc_p[4,] %>% select(n_gamma, p_fd, lb, ub, bb, non_bingers, bingers) %>% as.list()

names(args) <- c("mass", "p_fd", "lb", "ub", "bb", "non_bingers", "bingers")

args$mass <- args$mass[[1]]

args

rr_free_factories[5,]$current_fraction_factory[[1]](args)(250)




rr_clbr_factories <- makeCalibratedFactories(rr = r_cal, pc = pc_p, dh = dh_p)



fn <- calibrateSlope(target = .05, mass = function(x) dgamma(x, 1, 2), 0.03, 250)

fn



scen1[25,]$former_fraction[[1]](250)

scen1[925,]$current_fraction[[1]](c(125,250))





fn2 <- makeFormerCalibratedFactory()

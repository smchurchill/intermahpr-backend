library(ggplot2)
library(tidyverse)
library(magrittr)

RR <- rr_default
PC <- pc_default

EXT <- TRUE
lb <- LB <- 0.03
ub <- UB <- 250
BB <- c(53.8, 67.25)

IN_CUTS <- list("Female" = c(13.45, 26.9), "Male" = c(20.2, 40.4))

JOINED <- assemble(pc = PC, rr = RR, ext = EXT, lb = LB, ub = UB, bb = BB)
CUTTED <- compute_aafs(JOINED, IN_CUTS)

View(CUTTED)

JOINED %<>%
  mutate(CUTS = pmap(list(LB, IN_CUTS[GENDER], UB), c)) %>%
  mutate(CUMUL_F = map2(AAF_CMP, CUTS, ~.x(.y))) %>%
  mutate(AAF_GRP = map(CUMUL_F, ~diff(.x))) %>%
  mutate(AAF_TOT = map(CUMUL_F, ~`[[`(.x, length(.x))))

JOINED

JOINED[["AAF_CMP"]][[1]]

diff(
JOINED[["AAF_CMP"]][[1]](
JOINED[["CUTS"]][[1]])
)

for(level in names(CUTS)){
  print(level)
  print(CUTS[[level]])

}


RRF <- format_v0_rr(RR)
RRD <- derive_v0_rr(RRF, EXT)

PCF <- format_v0_pc(PC)

JOINT <- join_pc_rr(PCD, RRD)

JOINT

PIN <- 40

ggplot(data = data.frame(x=0), mapping = aes(x=x)) +
  stat_function(fun = function(x) c(diff(joint[[PIN, "AAF_CMP"]](x)), 0)) +
  xlim(0.03,250) +
  labs(title = paste(joint[[PIN, "OUTCOME"]],
                     "outcome due to",
                     joint[[PIN, "CONDITION"]],
                     joint[[PIN, "GENDER"]],
                     joint[[PIN, "AGE_GROUP"]],
                     joint[[PIN, "YEAR"]],
                     joint[[PIN, "REGION"]]),
       x = "average grams ethanol per day",
       y = "Cumulative AAF")

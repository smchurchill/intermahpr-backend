library(ggplot2)
library(tidyverse)
library(magrittr)

path = "C:\\Users\\samuelch\\Documents\\RStudio\\IMAHPROUTPUT\\"

dir <- tempdir()

dir

a <- intermahpr_raw()
b <- intermahpr_base(RelativeRisks = rr_default,
                     PrevalenceConsumption = pc_default,
                     FemaleLightModerateBarrier = 13.45,
                     FemaleModerateHeavyBarrier = 26.9,
                     FemaleBingeBarrier = 53.8,
                     MaleLightModerateBarrier = 20.2,
                     MaleModerateHeavyBarrier = 40.4,
                     MaleBingeBarrier = 67.25,
                     UpperBound = 250,
                     Extrapolation = TRUE,
                     OutputPath = path)

path

b

b[["InterMAHP_AAFs_morbidity"]]


dir(path)
dir("~/RStudio/IMAHPROUTPUT")

readr::write_csv(x = b[["InterMAHP_AAFs_morbidity"]], path = path)


check <- extract_prevcons(a)
check2 <- outcome_splitter(a, "Morbidity")

RR <- rr_default
PC <- pc_default

rrl <- levels(as.factor(RR$Gender))
pcl <- levels(as.factor(PC$Gender))


EXT <- TRUE
lb <- LB <- 0.03
ub <- UB <- 250
BB <- c(53.8, 67.25)

IN_CUTS <- list("Female" = c(13.45, 26.9), "Male" = c(20.2, 40.4))

cbl <- levels(as.factor(names(IN_CUTS)))

truths <- sapply(list(rrl, cbl, c("FEMALE", "MALE")), FUN = identical, pcl)

prod(truths)
sum(truths)

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

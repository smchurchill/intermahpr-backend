library(tidyverse)

RR <- rr_default
RRV <- format_v0_rr(RR)
RRV[, "EXT"] <- TRUE
RRV

RRB <- apply(RRV, 1, function(obs) compile_rr(obs))
RRB

RRV["CURVES"] <- RRB
RRV["CURVES"][[1]]

str(RRB)


PC <- pc_default
PCV <- format_v0_pc(PC)
PCS <- derive_v0_pc(PCV)
PCS[,-c(1,2,4:16)]



RRS[[1]][[1]][["GENDER"]]




D <- intermahpr::rr_default
E <- data.frame(D)
G <- as.tibble(E)
G

Fun <- function(row) {
  IM <- row[["IM"]]
  Gender <- row[["Gender"]]
  Function <- row[["Function"]]
  Outcome <- row[["Outcome"]]
  BingeF <- as.numeric(row[["BingeF"]])
  RR_FD <- as.numeric(row[["RR_FD"]])

  betas <- data.matrix(sapply(row[8:length(row)], as.numeric))
  intermahpr::simple_rr(betas,TRUE)
}

L <- deduce_relative_risk_curves_from_rr()

ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
  stat_function(fun = L[[54]][[4]]) +
  xlim(0.03, 150)

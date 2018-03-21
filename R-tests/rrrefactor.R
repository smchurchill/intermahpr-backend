library(tidyverse)

RR <- rr_default
RRV <- format_v0_rr(RR)
RRV[, "EXT"] <- TRUE
# RRV

RRV

allocate_as_fn <- function(in_data, var_name) {
  in_data[, var_name] <- list(list(var_name = function() var_name))
}

CURVE_LIST <- c("BASE_RR", "LNXT_RR", "BNGD_RR",
                "N_GAMMA", "INTGRND", "AAF_CMP")

RRR <- Reduce(f = allocate_as_fn, x = CURVE_LIST, init = RRV)

RRV[, "BASE_RR"] <- list(list("BASE_RR" = function(x) 0))
RRV[, "LNXT_RR"] <- list(list("LNXT_RR" = function(x) 0))
RRV[, "BNGD_RR"] <- list(list("BNGD_RR" = function(x) 0))

RRV[, c("BASE_RR", "LNXT_RR", "BNGD_RR")]

for(n in 1:nrow(RRV)) {
  RRV[n, "BASE_RR"] <- set_rr(RRV[n,])
  RRV[n, "CURVES"][["LNXT_RR"]] <- ext_rr(RRV[n,], set_rr(RRV[n,]))
  RRV[n, "CURVES"][["BNGD_RR"]] <- bng_rr(RRV[n,], set_rr(RRV[n,]))
}


RRB <- RRV %>%
  mutate()

  list(apply(RRV, 1, function(obs) compile_rr(obs)))
# RRB

RRV["CURVES"] <- RRB
RRV["CURVES"][[1]]

PC <- pc_default
PCV <- format_v0_pc(PC)
PCS <- derive_v0_pc(PCV)

JOINED <- PCS %>%
  full_join(RRV)

cmp <- compile_aaf(JOINED[1,])



PCS[,-c(1,2,4:16)]


for(n in 1:nrow(JOINED)) {
  DF <- JOINED[["DF"]][[n]]
  GAMMA_SHAPE <- JOINED[["GAMMA_SHAPE"]][[n]]
  GAMMA_SCALE <- JOINED[["GAMMA_SCALE"]][[n]]
  P_GAMMA <- function(x) {
    DF * dgamma(x, shape = GAMMA_SHAPE, scale = GAMMA_SCALE)
  }
  JOINED[["CURVES"]][[n]][["P_GAMMA"]] <- P_GAMMA
}

intgrnds <- apply(JOINED, 1, function(obs) int_factory(obs))

wng <- normalized_gamma_factory(JOINED[1,])




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

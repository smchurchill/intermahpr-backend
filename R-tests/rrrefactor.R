library(tidyverse)


RR <- rr_default
RRV <- format_v0_rr(RR)
RRV[, "EXT"] <- TRUE
# RRV

RRV



RRV <- add_column(RRV,
                  BASE_RR = zero,
                  LNXT_RR = zero,
                  BNGD_RR = zero)

for(n in 1:nrow(RRV)) {
  base <- set_rr(RRV[n,])
  lnxt <- ext_rr(RRV[n,], base)
  bngd <- bng_rr(RRV[n,], lnxt)
  RRV[[n, "BASE_RR"]] <- base
  RRV[[n, "LNXT_RR"]] <- lnxt
  RRV[[n, "BNGD_RR"]] <- bngd
}


PC <- pc_default
PCV <- format_v0_pc(PC)
PCS <- derive_v0_pc(PCV)

View(PCS)

JOINED <- PCS %>%
  full_join(RRV)

for(n in 1:nrow(JOINED)) {
  integrand <- integrand_factory(JOINED[n, ])
  JOINED[[n, "INTGRND"]] <- integrand
  aaf_computer <- aaf_factory(JOINED[n, ])
  JOINED[[n, "AAF_CMP"]] <- aaf_computer
  JOINED[[n, "AAF_FD"]] <- aaf_fd(JOINED[n, ])
}


View(JOINED)

PIN <- 200

ggplot(data = data.frame(x=0), mapping = aes(x=x)) +
  stat_function(fun = JOINED[[PIN, "INTGRND"]]) +
  xlim(0.03,250) +
  labs(title = paste(JOINED[[PIN, "CONDITION"]],
                     JOINED[[PIN, "GENDER"]],
                     JOINED[[PIN, "AGE_GROUP"]],
                     JOINED[[PIN, "YEAR"]],
                     JOINED[[PIN, "REGION"]]))

JOINED[1, "N_GAMMA"]

View(unnested)


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

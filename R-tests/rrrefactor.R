library(ggplot2)

RR <- rr_default
RRF <- format_v0_rr(RR)
RRD <- derive_v0_rr(RRF, TRUE)

PC <- pc_default
PCF <- format_v0_pc(PC)
PCD <- derive_v0_pc(PCF)

joined <- join_pc_rr(PCD[1,], RRD)



View(joined)

PCD

RRD[1:2,]

joint <- dplyr::full_join(PCD, RRD, by = "GENDER")

View(joint)

intgrnd <- intgrnd_factory(joint[1, ])
joint[[1, "INTGRND"]] <- intgrnd

aaf_cmp <- aaf_cmp_factory(joint[1, ])
joint[[1, "AAF_CMP"]] <- aaf_cmp



for(n in 1:nrow(joint)) {
  intgrnd <- intgrnd_factory(joint[n, ])
  joint[[n, "INTGRND"]] <- intgrnd
}

for(n in 1:nrow(joint)) {
  aaf_cmp <- aaf_cmp_factory(joint[n, ])
  joint[[n, "AAF_CMP"]] <- aaf_cmp
  #
  #     JOINT[[n, "AAF_FD"]] <- aaf_fd(JOINT[n, ])
}

joint

pts <- 1:100

joint1[[1, "AAF_CMP"]](pts)

JOINT <- join_pc_rr(PCD, RRD)

JOINT

PIN <- 95

ggplot(data = data.frame(x=0), mapping = aes(x=x)) +
  stat_function(fun = joint[[PIN, "AAF_CMP"]]) +
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

JOINED[1, "N_GAMMA"]

View(unnested)

integrand <- function(x) 3*(x^2)

pts <- 1:100

vector <- sapply(pts, function(x) integrate(f = integrand,
                                            lower = 0,
                                            upper = x)$value)


vector - (pts^3)


int <- function(x) {


}


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

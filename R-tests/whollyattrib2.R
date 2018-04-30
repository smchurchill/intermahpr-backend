library(tidyverse)
library(magrittr)



dh_in <- readr::read_csv(file.path("data-raw", "dh.csv"),
                         col_types = "??????dddd????????")


dh <- dh_in %>%
  rename(region = Province) %>%
  rename(condition = Condition_Alcohol)

dh$Outcome <- "Morbidity"

dh

DH <- dh %>%
  intermahpr::format_v0_dh() %>%
  filter(REGION == "10:BC" & YEAR == 2015)

DH$REGION = "BC"

DH

rr <- rr_default
RR <- rr %>%
  format_v0_rr() %>%
  derive_v0_rr(TRUE)

pc <- pc_default
PC <- pc_bc %>%
  format_v0_pc() %>%
  derive_v0_pc(
    bb = list("Female" = 50, "Male" = 60),
    lb = 0.03,
    ub = 250,
    gc = list("Female" = 1.258^2, "Male" = 1.171^2))

joined <- join_pc_rr(PC, RR)

joined

aaf <- base_aafs(joined)

aaf_split <- outcome_splitter(aaf)


SIMILAR <- c("IM", "REGION", "YEAR", "GENDER", "AGE_GROUP", "OUTCOME")

aaf_split$CONDITION <- NULL
aaf_split %<>% filter(OUTCOME == "Morbidity")
JOINED_FULL <- full_join(DH, aaf_split, by = SIMILAR)
JOINED_LEFT <- left_join(DH, aaf_split, by = SIMILAR)

JOINED_FULL %<>%
  mutate(BLOCK = as_factor(purrr::pmap_chr(
    list(REGION, YEAR, GENDER, AGE_GROUP, OUTCOME),
    .f = paste0))) %>%
  mutate(KEY = as_factor(purrr::pmap_chr(
    list(IM, BLOCK),
    .f = paste0)))

JCOPY <- JOINED_FULL

ims <- JCOPY$IM

ims[grepl("4.*[123]", ims)]
ims[grepl("(4.*[123]|6.*[15])", ims)]

ims[grepl("8.*5", ims)]


calibre <- filter(JOINED_FULL, grepl("(4.*[123]|6.*[15])", IM))

calibre$AAF_FD <- 0
calibre$AAF_TOTAL <- 1

for(n in 1:nrow(calibre)) {
  calibre[[n, "AAF_CMP"]] <- calibration_factory(calibre[n, ])
}


levels(JOINED_FULL$BLOCK)

for(block in levels(JOINED_FULL$BLOCK)) {
  curr_block <- filter(JOINED_FULL, BLOCK == block)

}

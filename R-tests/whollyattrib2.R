library(tidyverse)
library(magrittr)

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

aaf <- base_aafs(joined)

aaf_split <- outcome_splitter(aaf)

dh_in <- readr::read_csv(
  file.path("data-raw", "dh.csv"),
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

DH %<>% derive_v0_dh(PC)

DH

SIMILAR <- c("IM", "REGION", "YEAR", "GENDER", "AGE_GROUP", "OUTCOME")

aaf_split %<>%
  filter(OUTCOME == "Morbidity") %>%
  select(-c(CONDITION, DRINKERS, BB, LB, UB, N_GAMMA))

aaf_split

JOINED_FULL <- full_join(DH, aaf_split, by = SIMILAR)

JOINED_FULL %<>%
  mutate(BLOCK = as_factor(purrr::pmap_chr(
    list(REGION, YEAR, GENDER, AGE_GROUP, OUTCOME),
    .f = paste0))) %>%
  mutate(KEY = as_factor(purrr::pmap_chr(
    list(IM, REGION, YEAR, GENDER, AGE_GROUP, OUTCOME),
    .f = paste0))) %>%
  mutate(CC = substr(IM, 1, 3))

JCOPY <- JOINED_FULL

ims <- JCOPY$IM

ims[grepl("4.*[123]", ims)]
ims[grepl("(4.*[123]|6.*[15])", ims)]

ims[grepl("8.*5", ims)]


clbr <- filter(JOINED_FULL, grepl("(4.*[123]|5.*3|6.*[15])", IM))

clbr

clbr$AAF_FD <- 0
clbr$AAF_TOTAL <- 1

for(n in 1:nrow(calibre)) {
  calibre[[n, "AAF_CMP"]] <- calibration_factory(calibre[n, ])
}

check_calibre <- calibre %>%
  mutate(CUTS = pmap(list(0.03, 50, 100, 150, 200, 250), c)) %>%
  mutate(CUMUL_F = map2(AAF_CMP, CUTS, ~.x(.y))) %>%
  mutate(AAF_GRP = map(CUMUL_F, ~diff(.x))) %>%
  mutate(AAF_CD = unlist(map(CUMUL_F, ~`[`(.x, length(.x)))))

sims <- filter(JOINED_FULL, grepl("(4.*[467]|8...5|9...2)", IM))
sims_use <- filter(JOINED_FULL, grepl("(4.*5|all)", IM))

for(n in 1:nrow(sims)) {
  sims[[n, "AAF_CMP"]] <- scaled_aaf_cmp_factory(
    filter(sims_use, BLOCK == sims[[n, "BLOCK"]] & CC == sims[[n, "CC"]])
  )
}

check_sims <- sims %>%
  mutate(CUTS = pmap(list(0.03, 50, 100, 150, 200, 250), c)) %>%
  mutate(CUMUL_F = map2(AAF_CMP, CUTS, ~.x(.y))) %>%
  mutate(AAF_GRP = map(CUMUL_F, ~diff(.x))) %>%
  mutate(AAF_CD = unlist(map(CUMUL_F, ~`[`(.x, length(.x)))))

copy <- filter(JOINED_FULL, grepl("(5.*7|8...[12346]|9...[1345])", IM))
copy_use <- filter(JOINED_FULL, grepl("(6.*2|all)", IM))

for(n in 1:nrow(copy)) {
  category <- copy[[n, "CC"]]
  if(category == "(5)") {category = "(6)"}
  use <- filter(copy_use, BLOCK == copy[[n, "BLOCK"]] & CC == category)
  copy[[n, "AAF_CMP"]] <- copy_aaf_cmp_factory(use)
}

check_copy <- copy %>%
  mutate(CUTS = pmap(list(0.03, 50, 100, 150, 200, 250), c)) %>%
  mutate(CUMUL_F = map2(AAF_CMP, CUTS, ~.x(.y))) %>%
  mutate(AAF_GRP = map(CUMUL_F, ~diff(.x))) %>%
  mutate(AAF_CD = unlist(map(CUMUL_F, ~`[`(.x, length(.x)))))

KEEP <- filter(JOINED_FULL, !grepl("(4.*[123467]|5.*[37]|6.*[15]|[89].*\\()", IM))

COMBINED <- rbind(KEEP)

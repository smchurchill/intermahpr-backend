library(tidyverse)
library(magrittr)
library(pryr)

.data <- rr

names(.data)[grep("[0-9]$", names(.data))]


rr <- rr_default
RR <- rr %>%
  format_v1_rr() %>%
  derive_v1_rr(TRUE)

RR

pc <- pc_default
PC <- pc_bc %>%
  format_v1_pc %>%
  derive_v1_pc(
    bb = list("Female" = 50, "Male" = 60),
    lb = 0.03,
    ub = 250,
    gc = list("Female" = 1.258^2, "Male" = 1.171^2))

joined <- join_pc_rr(PC, RR)




aaf <- base_aafs(joined)

aaf_split <- outcome_splitter(joined)

dh_in <- readr::read_csv(
  file.path("data-raw", "dh.csv"),
  col_types = "??????dddd????????")


dh <- dh_in %>%
  rename(region = Province) %>%
  rename(condition = Condition_Alcohol) %>%
  filter(region == "10:BC" & Year == 2015)


dh$Outcome <- "Morbidity"

dh

DH <- dh %>%
  intermahpr::format_v1_dh() %>%
  filter(REGION == "10:BC" & YEAR == 2015)

DH$REGION = "BC"

DH %<>% derive_v1_dh(PC)

DH

res <- join_dh_aaf(DH, aaf_split)

res

check_res <- res %>%
  mutate(CUTS = pmap(list(0.03, 50, 100, 150, 200, 250), c)) %>%
  mutate(CUMUL_F = map2(AAF_CMP, CUTS, ~.x(.y))) %>%
  mutate(AAF_GRP = map(CUMUL_F, ~diff(.x))) %>%
  mutate(AAF_CD = map_dbl(CUMUL_F, ~`[`(.x, length(.x))))

check_res$AAF_CD

res[8,]$BLOCK

(res[8,]$ARGS)[[1]][["BLOCK"]]

str(res)










SIMILAR <- c("IM", "REGION", "YEAR", "GENDER", "AGE_GROUP", "OUTCOME")


aaf_split %<>%
  filter(OUTCOME == "Morbidity") %>%
  select(-c(CONDITION, DRINKERS, BB, LB, UB, N_GAMMA))

aaf_split

JOINED <- full_join(DH, aaf_split, by = SIMILAR)

JOINED %<>%
  mutate(BLOCK = as_factor(purrr::pmap_chr(
    list(REGION, YEAR, GENDER, AGE_GROUP, OUTCOME),
    .f = paste0))) %>%
  mutate(KEY = as_factor(purrr::pmap_chr(
    list(IM, REGION, YEAR, GENDER, AGE_GROUP, OUTCOME),
    .f = paste0))) %>%
  mutate(CC = substr(IM, 1, 3))

KEEP <- filter(JOINED,!grepl("(4...[^5]|5...[37]|6...[15]|[89]..\\()", IM))
MOD  <- filter(JOINED, grepl("(4...[^5]|5...[37]|6...[15]|[89]..\\([^a])", IM))

MOD

USE  <- filter(JOINED, grepl("(4.*5|6.*2|all)", IM))

c(address(MOD), refs(MOD))
tracemem(MOD)

for(n in 1:nrow(MOD)) {

  im <- MOD[[n, "IM"]]
  if(grepl("(4...[123]|5...3|6...[15])", im)) {
    expr <- NULL
    use <- MOD[n, ]
    get_aaf_fn <- calibration_factory
  }
  else {
    block <- MOD[[n, "BLOCK"]]
    cc <- MOD[[n, "CC"]]
    if(grepl("(4...[467]|8...5|9...2)", im)) {
      expr <- "(4.*5|all)"
      get_aaf_fn <- scaled_aaf_cmp_factory
    }
    else if(grepl("(5...7|8...[^5]|9...[^2])", im)) {
      if(grepl("5...7", im)) {
        cc <- "(6)"
      }
      expr <- "(6.*2|all)"
      get_aaf_fn <- copy_aaf_cmp_factory
    }
    else {
      stop(
        paste0(
          "IM of ",
          im,
          " does not match any of the following regular expressions:\n",
          "(4...[123]|5...3|6...[15])\n",
          "(4...[467]|8...5|9...2)\n",
          "(5...7|8...[^5]|9...[^2])",
          collapse = "\n"
        )
      )
    }
    use <- filter(USE, grepl(expr, IM) & BLOCK == block & CC == cc)
  }
  MOD[[n, "AAF_CMP"]] <- get_aaf_fn(use)
  MOD[[n, "AAF_FD"]] <- aaf_fd(use)
}

COMB <- rbind(KEEP, MOD)

check_comb <- COMB %>%
  mutate(CUTS = pmap(list(0.03, 50, 100, 150, 200, 250), c)) %>%
  mutate(CUMUL_F = map2(AAF_CMP, CUTS, ~.x(.y))) %>%
  mutate(AAF_GRP = map(CUMUL_F, ~diff(.x))) %>%
  mutate(AAF_CD = unlist(map(CUMUL_F, ~`[`(.x, length(.x))))) %>%
  mutate(AA_COUNTS = map2(AAF_GRP, COUNT, `*`))


clbr <- filter(JOINED_FULL, grepl("(4.*[123]|5.*3|6.*[15])", IM))

clbr

clbr$AAF_FD <- 0
clbr$AAF_TOTAL <- 1

for(n in 1:nrow(clbr)) {
  clbr[[n, "AAF_CMP"]] <- calibration_factory(clbr[n, ])
  clbr[[n, "AAF_FD"]] <- aaf_fd(clbr[n, ])
}

check_clbr <- clbr %>%
  mutate(CUTS = pmap(list(0.03, 50, 100, 150, 200, 250), c)) %>%
  mutate(CUMUL_F = map2(AAF_CMP, CUTS, ~.x(.y))) %>%
  mutate(AAF_GRP = map(CUMUL_F, ~diff(.x))) %>%
  mutate(AAF_CD = unlist(map(CUMUL_F, ~`[`(.x, length(.x)))))

sims <- filter(JOINED_FULL, grepl("(4.*[467]|8...5|9...2)", IM))
sims_use <- filter(JOINED_FULL, grepl("(4.*5|all)", IM))

for(n in 1:nrow(sims)) {
  use <- filter(sims_use, BLOCK == sims[[n, "BLOCK"]] & CC == sims[[n, "CC"]])
  sims[[n, "AAF_CMP"]] <- scaled_aaf_cmp_factory(use)
  sims[[n, "AAF_FD"]] <- aaf_fd(use)
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
  if(category == "(5)") {
    category = "(6)"
  }
  use <- filter(copy_use, BLOCK == copy[[n, "BLOCK"]] & CC == category)
  copy[[n, "AAF_CMP"]] <- copy_aaf_cmp_factory(use)
  copy[[n, "AAF_FD"]] <- aaf_fd(use)
}

check_copy <- copy %>%
  mutate(CUTS = pmap(list(0.03, 50, 100, 150, 200, 250), c)) %>%
  mutate(CUMUL_F = map2(AAF_CMP, CUTS, ~.x(.y))) %>%
  mutate(AAF_GRP = map(CUMUL_F, ~diff(.x))) %>%
  mutate(AAF_CD = unlist(map(CUMUL_F, ~`[`(.x, length(.x)))))

KEEP <- filter(JOINED_FULL, !grepl("(4.*[123467]|5.*[37]|6.*[15]|[89].*\\()", IM))

COMBINED <- rbind(KEEP)

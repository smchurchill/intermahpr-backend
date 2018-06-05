library(intermahpr)
library(tidyverse)
library(magrittr)

full_join(derive_pc(format_pc(pc_default)), derive_rr(format_rr(rr_default)))




rr_non_inj <- rr_default %>% filter(!grepl("all", IM))
rr_inj <- rr_default %>%
  filter(grepl("all", IM)) %>%
  mutate(CC = substr(IM, 1, 3)) %>%
  select(-IM)

dh <- readr::read_csv("data-raw/dh.csv")

conditions <- dh %>%
  select(IM, Condition_Alcohol) %>%
  distinct()

c_considered <- conditions %>%
  filter(!(grepl("A", IM)|grepl("NB-", Condition_Alcohol))) %>%
  mutate(IM = substr(IM, 1, 7)) %>%
  crossing(tibble(Gender = c("Male", "Female")))

non_inj_c <- c_considered %>%
  filter(!grepl("(8|9)", IM))
inj_c <- c_considered %>%
  filter(grepl("(8|9)", IM)) %>%
  mutate(CC = substr(IM, 1, 3)) %>%
  left_join(rr_inj) %>%
  select(-CC)




%>%
  left_join(rr)




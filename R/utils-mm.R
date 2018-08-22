#### Deaths/Hosps Specific Data Carpentry --------------------------------------

#' Prepare Death/Hospitalization Data
#'
#' @export
prepareMM <- function(.data) {
  message("Preparing morbidity and mortality input... ", appendLF = FALSE)

  .data %<>%
    clean(getExpectedVars("mm")) %>%
    collapseDeprecated()

  message("Done")

  .data
}


#' Collapse deprecated IM designations
#'
#'
collapseDeprecated <- function(.data) {
  .data %>%
    mutate(im = substring(text = im, first = 1, last = 7)) %>%
    group_by(region, year, gender, age_group, im, outcome) %>%
    summarise(count = sum(count)) %>%
    ungroup()
}

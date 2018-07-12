#' Format a scenario table for shiny output
#' @export
renderScenarioWide <- function(.data) {
  .data <- addFormerFraction(.data)
  .data <- addCurrentFraction(.data)
  .data <- addTotalFraction(.data)
  .data[c("current_fraction", "former_fraction")] <- NULL
  attr <-ifelse(
    .data$aaf == 0, 1,
    ((.data$attributability == "Wholly") / .data$aaf) +
      (.data$attributability == "Partially")
  )
  .data$aaf_cd <- .data$aaf_cd * attr
  .data$aaf <- .data$aaf * attr
  .data
}

renderScenarioLong <- function(.data) {


}

#' @export
runExample <- function() {
  appDir <- system.file("shiny-examples",
                        "intermahp-shiny",
                        package = "intermahpr")
  if(appDir == ""){
    stop("Could not find example directory.  Try reinstalling intermahpr.",
         call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}

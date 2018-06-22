#' Format a scenario table for shiny output
#' @export
formatForShinyOutput <- function(.data) {
  .data$aaf_fd <- computeFormerFraction(.data)
  .data$aaf_cd <- computeCurrentFraction(.data)
  .data$aaf <- computeTotalFraction(.data)
  .data[c("current_fraction", "former_fraction")] <- NULL
  attr <- ((.data$attributability == "Wholly") / .data$aaf) + (.data$attributability == "Partially")
  .data$aaf_cd <- .data$aaf_cd * attr
  .data$aaf <- .data$aaf * attr
  .data
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

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

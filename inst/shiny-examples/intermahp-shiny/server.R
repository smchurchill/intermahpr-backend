# Sam Churchill
# April 12 2018

# This is the server portion of a shiny demo for intermahp

server <- function(input, output, session) {
  shinyjs::hide("loadingContent")
  shinyjs::show("allContent")
}

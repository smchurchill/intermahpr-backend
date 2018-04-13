# Sam Churchill
# April 12 2018
#
# This is the ui portion of a shiny demo for intermahpr

function(request) {
  fluidPage(
    shinyjs::useShinyjs(),
    title = "InterMAHP",

    ## Add custom JS and CSS -

    ## Enclose the header in it's own section for nicer styling -
    div(
      id = "headerSection",
      h1("The International Model of Alcohol Harms and Policies"),

      ## Additional info -
      span(
        style = "font-size: 1.2em",
          span("CISUR package"),
        br(),

        span("April 12, 2018")
      )
    ),

    ## Initial Loading Screen -
    div(
      id = "loadingContent",
      h2("Loading...")
    ),

    ## Content goes here, hidden initially until full load -
    shinyjs::hidden(
      div(
        id = "allContent",

        ## sidebar - data input -
        sidebarLayout(
          sidebarPanel(
            h3("Relative Risk Data", style = "margin-top: 0;"),
            ## Upload your own RR csv or use a packaged one?
            selectInput(
              "source_rr",
              "",
              c(
                "Upload CSV" = "upload",
                "Use packaged relative risk curves" = "packaged"
              ),
              selected = "upload"
            ),
            ## Which IHD treatment among packaged risks? -
            conditionalPanel(
              "input.source_rr == 'packaged'",
              selectInput(
                "packaged_rr",
                "",
                c(
                  "Zhao IHD treatment" = "zhao_rr",
                  "Roerecke IHD treatment" = "roerecke_rr"
                ),
                selected = "zhao_rr"
              )
            ),
            ## File upload for RR -
            conditionalPanel(
              "input.source_rr == 'upload'",
              fileInput(
                "uploaded_rr",
                "",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values",
                  "text/plain",
                  ".csv"
                ),
                buttonLabel = "Choose File",
                placeholder = "RelRisks.csv"
              )
            ),
            br(),

            h3("Prevalence and Consumption Data", style = "margin-top: 0;"),
            ## Upload your own PC csv or use a packaged one? -
            selectInput(
              "source_pc",
              "",
              c(
                "Upload CSV" = "upload",
                "Use packaged prevalence & consumption" = "packaged"
              ),
              selected = "upload"
            ),
            ## Which PC data among packaged? -
            conditionalPanel(
              "input.source_pc == 'packaged'",
              selectInput(
                "packaged_pc",
                "",
                c(
                  "British Columbia, 2015" = "pc_bc",
                  "Canada, 2015" = "pc_can",
                  "Combined" = "pc_bc_can"
                ),
                selected = "pc_bc"
              )
            ),
            ## File upload for PC -
            conditionalPanel(
              "input.source_pc == 'upload'",
              fileInput(
                "uploaded_pc",
                "",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values",
                  "text/plain",
                  ".csv"
                ),
                buttonLabel = "Choose File",
                placeholder = "PrevCons.csv"
              )
            ),
            h3("Drinking Group Definitions"),
            p("All measurements are made in units of grams-ethanol per day."),
            ## Tabs for numerical drinking group defintions -
            tabsetPanel(
              id = "numericalParamTab",
              type = "tabs",

              ## Tab accepting global parameter input (methods) -
              tabPanel(
                title = "Methods",
                id = "methodTab",
                numericInput(
                  "upper_bound",
                  "Upper limit of consumption",
                  value = 250,
                  min = 0,
                  step = 1
                ),
                br(),
                selectInput(
                  "ext",
                  "Dose response extrapolation method",
                  c(
                    "Linear" = TRUE,
                    "Capped" = FALSE
                  ),
                  selected = TRUE
                )
              ),
              ## Tab accepting female drinking group definitions -
              tabPanel(
                title = "Female",
                id = "femaleTab",
                numericInput(
                  "lm_f",
                  "Light-moderate barrier",
                  value = 15,
                  min = 0,
                  step = 1
                ),
                numericInput(
                  "mh_f",
                  "moderate-heavy barrier",
                  value = 30,
                  min = 0,
                  step = 1
                ),
                numericInput(
                  "bb_f",
                  "Binge barrier",
                  value = 50,
                  min = 0,
                  step = 1
                )
              ),
              ## Tab accepting male drinking group definitions -
              tabPanel(
                title = "Male",
                id = "maleTab",
                numericInput(
                  "lm_m",
                  "Light-moderate barrier",
                  value = 20,
                  min = 0,
                  step = 1
                ),
                numericInput(
                  "mh_m",
                  "Moderate-heavy barrier",
                  value = 40,
                  min = 0,
                  step = 1
                ),
                numericInput(
                  "bb_m",
                  "Binge barrier",
                  value = 65,
                  min = 0,
                  step = 1
                )
              )
            ) ## end tabsetPanel -
          ), ## end sidebarPanel -
          mainPanel(
            wellPanel(
              tabsetPanel(
                id = "resultsTab",
                type = "tabs",

                tabPanel(
                  title = "Show data",
                  id = "tableTab"
                )
                ## To be added -
              )
            )
          )
        )
      )
    ) ## end shinyjs::hidden -
  ) ## end fluidPage -
}

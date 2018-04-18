# Sam Churchill
# April 12 2018
#
# This is the ui portion of a shiny demo for intermahpr

function(request) {
  fluidPage(
    shinyjs::useShinyjs(),
    title = "InterMAHP",

    ## Add custom JS and CSS ----
    shiny::singleton(
      tags$head(
        includeCSS(file.path("www", "intermahp.css"))
      )
    ),

    ## Enclose the header in it's own section ----
    ## for nicer styling
    div(
      id = "headerSection",
      h1("The International Model of Alcohol Harms and Policies"),

      ## Additional info ----
      span(
        style = "font-size: 1.2em",
          span("CISUR package"),
        br(),

        span("April 12, 2018")
      )
    ),

    ## Initial Loading Screen ----
    div(
      id = "loadingContent",
      h2("Loading...")
    ),

    ## Content goes here, hidden initially until full load ----
    shinyjs::hidden(
      div(
        id = "allContent",

        ## sidebar - data input ----
        sidebarLayout(
          sidebarPanel(
            h3("Relative Risk Data", style = "margin-top: 0;"),
            ## Upload your own RR csv or use a packaged one? ----
            selectInput(
              "source_rr",
              "",
              c(
                "Upload CSV" = "upload",
                "Use packaged" = "packaged"
              ),
              selected = "packaged"
            ),
            ## Which IHD treatment among packaged risks? ----
            conditionalPanel(
              style = "height: 94px;",
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
            ## File upload for RR ----
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
            ## Upload your own PC csv or use a packaged one? ----
            selectInput(
              "source_pc",
              "",
              c(
                "Upload CSV" = "upload",
                "Use packaged" = "packaged"
              ),
              selected = "packaged"
            ),
            ## Which PC data among packaged? ----
            conditionalPanel(
              style = "height: 74px;",
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
            ## File upload for PC ----
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
            p("All measurements are in units of grams-ethanol per day."),
            ## Tabs for numerical drinking group defintions ----
            tabsetPanel(
              id = "paramTabset",
              type = "tabs",

              ## Tab accepting global parameter input (methods) ----
              tabPanel(
                title = "Methods",
                id = "methodTab",
                br(),
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
              ## Tab accepting female drinking group definitions ----
              tabPanel(
                title = "Female",
                id = "femaleTab",
                br(),
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
              ## Tab accepting male drinking group definitions ----
              tabPanel(
                title = "Male",
                id = "maleTab",
                br(),
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
            ) ## end tabsetPanel ----
          ), ## end sidebarPanel ----
          ## Main Panel ----
          mainPanel(
            tabsetPanel(
              id = "resultsTabset",
              type = "tabs",

              ## Groups Display Tab ----
              tabPanel(
                title = "Drinking Groups",
                id = "drinkgrpTab",
                br(),
                fluidRow(
                  ## Female Grouping Display ----
                  column(
                    width = 6,
                    class = "scrollable",
                    h4("Female"),
                    tableOutput("drinking_groups_female"),
                    textOutput("binge_barrier_female")
                  ),
                  ## Male Grouping Display ----
                  column(
                    width = 6,
                    class = "scrollable",
                    h4("Male"),
                    tableOutput("drinking_groups_male"),
                    textOutput("binge_barrier_male")
                  )
                )
              ) ## end Groups display ----
              ,
              ## Relative Risk Display ----
              tabPanel(
                title = "Relative Risk",
                id = "relriskTab",
                br(),
                tabsetPanel(
                  id = "relrisksTabset",
                  type = "tabs",
                  tabPanel(
                    title = "Table",
                    id = "relrisksTable",
                    radioButtons(
                      "rr_table_type",
                      "",
                      choices = c("Raw" = "raw", "Formatted" = "for"),
                      selected = "raw",
                      inline = TRUE
                    ),
                    DT::dataTableOutput("rrTable")
                  ),
                  tabPanel(
                    title = "Plots",
                    id = "relrisksPlots",
                    fluidRow(
                      br(),
                      column(
                        width = 4,
                        uiOutput("select_condition_rr")
                      ),
                      column(
                        width = 4,
                        uiOutput("select_gender_rr")
                      ),
                      column(
                        width = 4,
                        uiOutput("select_outcome_rr")
                      )
                    )
                  )
                )
              ) ## end Rel Risks display ----
              ,
              ## PrevCons display ----
              tabPanel(
                title = "Prevalence and Consumption",
                id = "prevconsTab",
                br(),
                tabsetPanel(
                  id = "prevconsTabset",
                  type = "tabs",
                  tabPanel(
                    title = "Table",
                    id = "prevconsTable",
                    radioButtons(
                      "inprev_display_type",
                      "",
                      choices = c("Decimal" = "dec", "Percentage" = "per"),
                      selected = "per",
                      inline = TRUE
                    ),
                    DT::dataTableOutput("pcTable")
                  ),
                  tabPanel(
                    title = "Plots",
                    id = "prevconsPlots",
                    fluidRow(
                      br(),
                      column(
                        width = 6,
                        uiOutput("select_gender_pc")
                      ),
                      column(
                        width = 6,
                        uiOutput("select_age_group_pc")
                      )
                    )
                  )
                )
              ) ## end PrevCons display ----
              ,
              ## AAF display ----
              tabPanel(
                title = "Alcohol Attributable Fractions",
                id = "aafTab",
                br(),
                tabsetPanel(
                  id = "aafTabset",
                  type = "tabs",
                  tabPanel(
                    title = "Prevalence and Consumption",
                    id = "aafPCTable",
                    radioButtons(
                      "outprev_display_type",
                      "",
                      choices = c("Decimal" = "dec", "Percent" = "per"),
                      selected = "per",
                      inline = TRUE
                    ),
                    DT::dataTableOutput("prev_cons_output")
                  ),
                  tabPanel(
                    title = "Mortality",
                    id = "mortTable",
                    radioButtons(
                      "mort_aaf_display_type",
                      "",
                      choices = c("Decimal" = "dec", "Percent" = "per"),
                      selected = "per",
                      inline = TRUE
                    ),
                    DT::dataTableOutput("mortality_aaf")
                  ),
                  tabPanel(
                    title = "Morbidity",
                    id = "morbTable",
                    radioButtons(
                      "morb_aaf_display_type",
                      "",
                      choices = c("Decimal" = "dec", "Percent" = "per"),
                      selected = "per",
                      inline = TRUE
                    ),
                    DT::dataTableOutput("morbidity_aaf")
                  ),
                  tabPanel(
                    title = "Plots",
                    id = "aafPlots"
                  )
                )
              )## end AAF Display ----
              ,
              ## Harm Reduction Display ----
              tabPanel(
                title = "Harm Reduction",
                id = "harmReduction"
              )## end Harm Reduction Display ----
            )
          )## end mainPanel ----
        )
      )
    ) ## end shinyjs::hidden ----
  ) ## end fluidPage ----
}

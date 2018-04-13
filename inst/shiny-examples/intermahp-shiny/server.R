# Sam Churchill
# April 12 2018

# This is the server portion of a shiny demo for intermahp

server <- function(input, output, session) {

  ## Hide/show content ----
  shinyjs::hide("loadingContent")
  shinyjs::show("allContent")

  ## Female drinking groups ----
  output$drinking_groups_female <- renderTable({
    tibble::data_frame(Group = c("Light", "Moderate", "Heavy"),
               From = c(0.03, input$lm_f, input$mh_f),
               To = as.double(c(input$lm_f, input$mh_f, input$upper_bound))
    )
  })

  output$binge_barrier_female <- renderText({
    paste("Binge level: at least",
          input$bb_f,
          "grams-ethanol in a single drinking event."
    )
  })

  ## Male drinking groups ----
  output$drinking_groups_male <- renderTable({
    tibble::data_frame(Group = c("Light", "Moderate", "Heavy"),
               From = c(0.03, input$lm_m, input$mh_m),
               To = as.double(c(input$lm_m, input$mh_m, input$upper_bound))
    )
  })

  output$binge_barrier_male <- renderText({
    paste("Binge level: at least",
          input$bb_m,
          "grams-ethanol in a single drinking event."
    )
  })

  ## Relative risk data input ----
  rrData <- reactive({
    ds <- NULL
    if(input$source_rr == "packaged") {
      ds <- switch(input$packaged_rr,
                   "zhao_rr" = intermahpr::rr_zhao,
                   "roerecke_rr" = intermahpr::rr_roerecke
      )
    }
    else {
      inFile <- input$uploaded_rr
      if(!is.null(inFile)){
        ds <- read.csv(inFile$datapath)
      }
    }
    ds
  })

  frrData <- reactive({
    ds <- NULL
    if(!is.null(rrData())) {
      ds <- format_v0_rr(rrData())
    }
    ds
  })

  drrData <- reactive({
    ds <- NULL
    if(!is.null(frrData())) {
      ds <- derive_v0_rr(rr = frrData(), ext = input$ext)
    }
    ds
  })

  ## Prev cons data input ----
  pcData <- reactive({
    ds <- NULL
    if(input$pc_type == "packaged") {
      ds <- switch(input$packaged_pc,
                   "pc_bc" = intermahpr::pc_bc,
                   "pc_can" = intermahpr::pc_can,
                   "pc_bc_can" = intermahpr::pc_bc_can
      )
    }
    else {
      inFile <- input$uploaded_pc
      if(!is.null(inFile)){
        ds <- read.csv(inFile$datapath)
      }
    }
    ds
  })

  fpcData <- reactive({
    ds <- NULL
    if(!is.null(pcData())) {
      ds <- format_v0_pc(pcData())
    }
    ds
  })

  dpcData <- reactive({
    ds <- NULL
    if(!is.null(fpcData())) {
      ds <- derive_v0_pc(pc = fpcData(),
                         bb = list("Female" = input$bb_f, "Male" = input$bb_m),
                         lb = 0.03,
                         ub = input$upper_bound,
                         gc = list("Female" = 1.258^2, "Male" = 1.171^2))
    }
  })

  output$pcTable <- renderDataTable(
    {
      pcData()
    },
    options = list(pageLength = 18,
                   lengthMenu = c(12, 18, 36, 72))
  )

  ## Relative risk table render ----
  output$rrTable <- DT::renderDataTable(
    {
      ds <- NULL
      if(input$rr_table_type == "raw") {
        ds <- rrData()
      }
      else {
        ds <- frrData()
      }
      ds
    },
    filter = 'top',
    server = FALSE,
    options = list(pageLength = 10,
                   scrollX = TRUE,
                   autoWidth = TRUE
    )
  )

  ## Relative risk plot render ----
  conditions <- reactive({
    dplyr::distinct(frrData()["CONDITION"])$CONDITION
  })

  genders_rr <- reactive({
    distinct(frrData$GENDER)
  })

  outcomes <- reactive({
    distinct(frrData$OUTCOME)
  })

  output$select_condition <- renderUI({
    selectInput(
      "condition_rr_plot",
      "Conditions",
      choices = conditions(),
      multiple = TRUE,
      selectize = TRUE
    )
  })

  output$select_gender <- renderUI({
  })

  output$select_outcome <- renderUI({
  })

  ## Prev Cons table render ----

}

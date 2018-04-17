# Sam Churchill
# April 12 2018

# This is the server portion of a shiny demo for intermahp

server <- function(input, output, session) {
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
    if(input$source_pc == "packaged") {
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

  ## Relative risk table render ----
  output$rrTable <- DT::renderDataTable(
    {
      if(is.null(rrData())){return(NULL)}

      if(input$rr_table_type == "raw") {
        ds <- rrData()
      }
      else {
        ds <- frrData()
      }
      ds <- DT::datatable(
        data = ds,
        options = list(
          pageLength = 12,
          lengthMenu = c(12,18,36,72),
          scrollX = TRUE,
          autoWidth = TRUE,
          columnDefs = list(
            list(
              targets = 0,
              visible = FALSE
            ),
            list(
              targets = 2,
              render = DT::JS(
                "function(data, type, row, meta) {",
                "return type === 'display' && data.length > 15 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 13) + ",
                "'...</span>' : data;",
                "}")
            )
          ),
          filter = 'top'
        )
      )
      ds
    }
  )

  ## Relative risk plot render ----
  conditions_rr <- reactive({
    dplyr::distinct(frrData()["CONDITION"])$CONDITION
  })

  genders_rr <- reactive({
    dplyr::distinct(frrData()["GENDER"])$GENDER
  })

  outcomes_rr <- reactive({
    dplyr::distinct(frrData()["OUTCOME"])$OUTCOME
  })

  output$select_condition_rr <- renderUI({
    selectInput(
      "condition_rr_plot",
      "Conditions",
      choices = conditions_rr(),
      multiple = TRUE,
      selectize = TRUE
    )
  })

  output$select_gender_rr <- renderUI({
  })

  output$select_outcome_rr <- renderUI({
  })

  ## Prev Cons table render ----
  output$pcTable <- DT::renderDataTable(
    {
      if(is.null(pcData())) {return(NULL)}
      if(input$pc_table_type == "raw") {
        table <- pcData()
      }
      else{
        table <- fpcData()
      }
      DT::datatable(
        data = table,
        options = list(
          filter = 'top',
          pageLength = 18,
          lengthMenu = c(12,18,36,72),
          scrollX = TRUE,
          columnDefs = list(
            list(
              targets = 0,
              visible = FALSE
            )
          )
        )
      )
    }
  )

  ## Prev Cons plot render ----
  genders_pc <- reactive({
    dplyr::distinct(fpcData()["GENDER"])$GENDER
  })

  output$select_gender_pc <- renderUI({
    selectInput(
      "genders_pc_plot",
      "Genders",
      choices = genders_pc(),
      multiple = TRUE,
      selectize = TRUE
    )
  })

  output$select_gender_rr <- renderUI({
  })

  output$select_outcome_rr <- renderUI({
  })

  ## AAF objects ----
  assembled_fn <- reactive({
    ds <- NULL
    if((!is.null(dpcData())) && (!is.null(drrData()))) {
      ds <- join_pc_rr(pc = dpcData(), rr = drrData())
    }
    ds
  })

  cuts <- reactive({
    list("Female" = c(input$lm_f, input$mh_f),
         "Male" = c(input$lm_m, input$mh_m))
  })

  ready_for_evaluation <- reactive({
    ds <- NULL
    if(!is.null(assembled_fn())) {
      ds <- add_cutpoints(aaf_table = assembled_fn(), cuts = cuts())
    }
    ds
  })

  prev_cons_output_table <- reactive({
    ds <- NULL
    if(!is.null(ready_for_evaluation())) {
      ds <- intermahpr::extract_prevcons(aaf_table = ready_for_evaluation())
    }
    ds
  })

  output$prev_cons_output <- DT::renderDataTable(
    {
      DT::datatable(
        data = prev_cons_output_table(),
        options = list(
          filter = 'top',
          pageLength = 18,
          lengthMenu = c(12,18,36,72),
          scrollX = TRUE,
          columnDefs = list(
            list(
              targets = 0,
              visible = FALSE
            )
          )
        )
      )
    }
  )

  evaluated <- reactive({
    ds <- NULL
    if(!is.null(ready_for_evaluation())) {
      ds <- evaluate_at_cutpoints(ready_for_evaluation())
    }
    ds
  })

  mort_aaf_table <- reactive({
    ds <- NULL
    if(!is.null(evaluated())) {
      ds <- outcome_splitter(aaf_table = evaluated(), outcome = "Mortality")
    }
    ds
  })

  morb_aaf_table <- reactive({
    ds <- NULL
    if(!is.null(evaluated())) {
      ds <- outcome_splitter(aaf_table = evaluated(), outcome = "Morbidity")
    }
    ds
  })

  output$mortality_aaf <- DT::renderDataTable(
    {
      if(is.null(mort_aaf_table())) { return(NULL) }
      table <- DT::datatable(
        data = mort_aaf_table(),
        options = list(
          pageLength = 12,
          lengthMenu = c(12,18,36,72),
          scrollX = TRUE,
          autoWidth = TRUE,
          columnDefs = list(
            list(
              targets = 0,
              visible = FALSE
            ),
            list(
              targets = 6,
              render = DT::JS(
                "function(data, type, row, meta) {",
                  "return type === 'display' && data.length > 15 ?",
                  "'<span title=\"' + data + '\">' + data.substr(0, 13) + ",
                  "'...</span>' : data;",
                "}")
              )
            ),
          filter = 'top'
        )
      )
      if(input$mort_aaf_display_type == "dec") {
        table <- DT::formatRound(table = table,
                                 columns = c("AAF_FD", "AAF_LD", "AAF_MD",
                                             "AAF_HD", "AAF_TOTAL"),
                                 digits = 4)
      }
      else {
        table <- DT::formatPercentage(table = table,
                                      columns = c("AAF_FD", "AAF_LD", "AAF_MD",
                                             "AAF_HD", "AAF_TOTAL"),
                                      digits = 2)
      }
      table
    }
  )

  output$morbidity_aaf <- DT::renderDataTable(
    {
      if(is.null(morb_aaf_table())) { return(NULL) }
      table <- DT::datatable(
        data = morb_aaf_table(),
        options = list(
          pageLength = 12,
          lengthMenu = c(12,18,36,72),
          scrollX = TRUE,
          autoWidth = TRUE,
          columnDefs = list(
            list(
              targets = 0,
              visible = FALSE
            ),
            list(
              targets = 6,
              render = DT::JS(
                "function(data, type, row, meta) {",
                "return type === 'display' && data.length > 15 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 13) + ",
                "'...</span>' : data;",
                "}")
            )
          ),
          filter = 'top'
        )
      )
      if(input$morb_aaf_display_type == "dec") {
        table <- DT::formatRound(table = table,
                                 columns = c("AAF_FD", "AAF_LD", "AAF_MD",
                                             "AAF_HD", "AAF_TOTAL"),
                                 digits = 4)
      }
      else {
        table <- DT::formatPercentage(table = table,
                                      columns = c("AAF_FD", "AAF_LD", "AAF_MD",
                                                  "AAF_HD", "AAF_TOTAL"),
                                      digits = 2)
      }
      table
    }
  )


  ## Hide/show content ----
  shinyjs::hide("loadingContent")
  shinyjs::show("allContent")
}

# Sam Churchill
# April 12 2018
#
# This is the server portion of a shiny demo for intermahp

source("global.R") # import functions
source("helper.R") # Have helpers available


#'
#'@importFrom dplyr filter
#'@importFrom magrittr %>% %<>%
#'

server <- function(input, output, session) {
  ## reactive() Objects ----

  ## * Relative Risk ----
  rrData <- reactive({
    if(input$source_rr == "upload"){
      inFile <- input$uploaded_rr
      if(is.null(inFile)) {return(NULL)}
      ds <- read.csv(inFile$datapath)
    } else {
      ds <- switch(input$packaged_rr,
                   "zhao_rr" = intermahpr::rr_zhao,
                   "roerecke_rr" = intermahpr::rr_roerecke
      )
    }
    ds
  })

  frrData <- reactive(format_v0_rr(rrData()))
  drrData <- reactive(derive_v0_rr(rr = frrData(), ext = input$ext))

  conditions_rr <- reactive({dplyr::distinct(frrData()["CONDITION"])$CONDITION})
  genders_rr <- reactive({dplyr::distinct(frrData()["GENDER"])$GENDER})
  outcomes_rr <- reactive({dplyr::distinct(frrData()["OUTCOME"])$OUTCOME})


  ## * Prevalence and Consumption ----
  pcData <- reactive({
    if(input$source_pc == "upload") {
      inFile <- input$uploaded_pc
      if(is.null(inFile)) {return(NULL)}
      ds <- read.csv(inFile$datapath)
    } else {
      ds <- switch(input$packaged_pc,
                   "pc_bc" = intermahpr::pc_bc,
                   "pc_can" = intermahpr::pc_can,
                   "pc_bc_can" = intermahpr::pc_bc_can
      )
    }
    ds
  })

  fpcData <- reactive(format_v0_pc(pcData()))
  dpcData <- reactive(
    derive_v0_pc(
      pc = fpcData(),
      bb = list("Female" = input$bb_f, "Male" = input$bb_m),
      lb = 0.03,
      ub = input$upper_bound,
      gc = list("Female" = 1.258^2, "Male" = 1.171^2)
    )
  )

  regions_pc <- reactive({dplyr::distinct(fpcData()["REGION"])$REGION})
  years_pc <- reactive({dplyr::distinct(fpcData()["YEAR"])$YEAR})
  genders_pc <- reactive({dplyr::distinct(fpcData()["GENDER"])$GENDER})
  agegroups_pc <- reactive({dplyr::distinct(fpcData()["AGE_GROUP"])$AGE_GROUP})

  ## * AAF ----
  assembled_fn <- reactive({
    if(is.null(dpcData()) || is.null(drrData())) {return(NULL)}
    join_pc_rr(pc = dpcData(), rr = drrData())
  })

  cuts <- reactive({
    list("Female" = c(input$lm_f, input$mh_f),
         "Male" = c(input$lm_m, input$mh_m))
  })

  ready_for_evaluation <- reactive({
    if(is.null(assembled_fn())) {return(NULL)}
    add_cutpoints(aaf_table = assembled_fn(), cuts = cuts())
  })

  prev_cons_output_table <- reactive({
    if(is.null(ready_for_evaluation())) {return(NULL)}
    extract_prevcons(aaf_table = ready_for_evaluation())
  })

  evaluated <- reactive({
    if(is.null(ready_for_evaluation())) {return(NULL)}
    evaluate_at_cutpoints(ready_for_evaluation())
  })

  mort_aaf_table <- reactive({
    if(is.null(evaluated())) {return(NULL)}
    outcome_splitter(aaf_table = evaluated(), outcome = "Mortality")
  })

  morb_aaf_table <- reactive({
    if(is.null(evaluated())) {return(NULL)}
    outcome_splitter(aaf_table = evaluated(), outcome = "Morbidity")
  })

  ## Output ----

  ## * datatables ----
  ## ** arguments ----
  output$rrTable <- DT::renderDataTable(
    {
      if(is.null(rrData())) {return(NULL)}

      options <- base_options
      options$columnDefs <- list(cd_hide0,
                                 cd_dots(2, 15, 13))

      if(input$rr_table_type == "raw") {
        ds <- rrData()
      }
      else {
        ds <- frrData()
      }

      DT::datatable(
        data = ds,
        filter = "top",
        extensions = "Buttons",
        options = options
      )
    }
  )

  output$pcTable <- DT::renderDataTable(
    {
      if(is.null(pcData())) {return(NULL)}

      options <- base_options
      options$columnDefs <- list(cd_hide0)

      table <- DT::datatable(
        data = fpcData(),
        filter = "top",
        extensions = "Buttons",
        options = options
      )

      formatColumns <- c("P_LA", "P_FD", "P_CD", "P_BD")

      reformat(tbl = table,
               type = input$inprev_display_type,
               col = formatColumns)
    }
  )

  ## ** results ----
  output$prev_cons_output <- DT::renderDataTable(
    {
      if(is.null(pcData())) {return(NULL)}

      options <- base_options
      options$columnDefs <- list(cd_hide0)

      table <- DT::datatable(
        data = prev_cons_output_table(),
        filter = "top",
        extensions = "Buttons",
        options = options
      )

      table <- DT::formatRound(
        table = table,
        columns = c("PCC_AMONG_DRINKERS", "GAMMA_SHAPE", "GAMMA_SCALE"),
        digits = 3
      )

      formatColumns <- c("P_LA", "P_FD", "P_CD", "P_LD", "P_MD", "P_HD")

      reformat(tbl = table,
               type = input$outprev_display_type,
               col = formatColumns)
      }
  )

  output$mortality_aaf <- DT::renderDataTable(
    {
      if(is.null(mort_aaf_table())) { return(NULL) }

      options <- base_options
      options$columnDefs <- list(cd_hide0,
                                 cd_dots(6, 15, 13))

      table <- DT::datatable(
        data = mort_aaf_table(),
        filter = "top",
        extensions = "Buttons",
        options = options
      )

      formatColumns <- c("AAF_FD", "AAF_LD", "AAF_MD", "AAF_HD", "AAF_TOTAL")

      reformat(tbl = table,
               type = input$mort_aaf_display_type,
               col = formatColumns)
    }
  )

  output$morbidity_aaf <- DT::renderDataTable(
    {
      if(is.null(morb_aaf_table())) { return(NULL) }

      options <- base_options
      options$columnDefs <- list(cd_hide0,
                                 cd_dots(6, 15, 13))

      table <- DT::datatable(
        data = morb_aaf_table(),
        filter = "top",
        extensions = "Buttons",
        options = options
      )

      formatColumns <- c("AAF_FD", "AAF_LD", "AAF_MD", "AAF_HD", "AAF_TOTAL")

      reformat(tbl = table,
               type = input$morb_aaf_display_type,
               col = formatColumns)
    }
  )

  ## * tables ----
  ## ** Female drinking groups ----
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

  ## ** Male drinking groups ----
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

  ## * renderUI ----

  ## ** RR Plot selectors ----
  output$select_condition_rr <- renderUI({
    selectInput(
      "conditions_to_plot_rr",
      "Conditions",
      choices = conditions_rr(),
      multiple = TRUE,
      selectize = TRUE
    )
  })

  output$select_gender_rr <- renderUI({
    selectInput(
      "genders_to_plot_rr",
      "Genders",
      choices = genders_rr(),
      multiple = TRUE,
      selectize = TRUE
    )
  })

  output$select_outcome_rr <- renderUI({
    selectInput(
      "outcomes_to_plot_rr",
      "Outcomes",
      choices = outcomes_rr(),
      multiple = TRUE,
      selectize = TRUE
    )
  })

  ## ** PC Plot selectors ----
  output$select_region_pc <- renderUI({
    selectInput(
      "regions_to_plot_pc",
      "Regions",
      choices = regions_pc(),
      multiple = TRUE,
      selectize = TRUE
    )
  })

  output$select_year_pc <- renderUI({
    selectInput(
      "years_to_plot_pc",
      "Years",
      choices = years_pc(),
      multiple = TRUE,
      selectize = TRUE
    )
  })

  output$select_gender_pc <- renderUI({
    selectInput(
      "genders_to_plot_pc",
      "Genders",
      choices = genders_pc(),
      multiple = TRUE,
      selectize = TRUE
    )
  })

  output$select_agegroups_pc <- renderUI({
    selectInput(
      "agegroups_to_plot_pc",
      "Age Groups",
      choices = agegroups_pc(),
      multiple = TRUE,
      selectize = TRUE
    )
  })

  ## * Plots ----

  ## ** RR ----
  rrPlotPts <- reactive({
    input$updateRRPlotList
    ds <- rrPlotList()
    evalAt <- seq_len(500) * input$upper_bound / 500
    ds
    lval <- dplyr::mutate(lval = map(rrPlotList(), evalAt, ~.x(evalAt)))
    bval <- dplyr::mutate()
  })

  rrPlotList <- reactive({
    ds <- drrData()
    ds %<>%
      dplyr::filter(
        CONDITION %in% input$conditions_to_plot_rr &
          GENDER %in% input$genders_to_plot_rr &
          OUTCOME %in% input$outcomes_to_plot_rr
      )
    ds[,1:7]
  })

  output$test <- renderTable(rrPlotList())

  ## Hide/show content ----
  shinyjs::hide("loadingContent")
  shinyjs::show("allContent")
}

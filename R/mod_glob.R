#' glob UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_glob_ui <- function(id) {
  ns <- NS(id)

  tabsetPanel(
    id = ns("glob_tabs"),
    tabPanel("DAList creation",
             value = "glob_dalist_creation",
             mod_glob_DAlist_ui(ns("DAlist")),
             ### Create an actionButton at the end of DAlist creation step to unlock next step
             fluidRow(
               column(12, offset = 0, style = "text-align: right;",
                      h4("Step 5: Move on the next step:"),
                      actionButton(ns("DAlist_done"), "GO to the Filter sample(s) step",
                                   class = "btn-lg btn-success")
               )
             )),

    tabPanel("Filter sample(s)",
             value = "glob_filter_sample",
             mod_glob_filt_sample_ui(ns("filt_sample")),
             ### Create an actionButton at the end of filter_sample(s) step to unlock next step
             fluidRow(
               column(12, offset = 0, style = "text-align: right;",
                      h4("Step 4: Move on the next step:"),
                      actionButton(ns("filter_sample_done"), "GO to the Filter protein(s) step",
                                   class = "btn-lg btn-success")
               )
             )),

    tabPanel("Filter proteins",
             value = "glob_filter_protein",
             mod_glob_filt_prot_ui(ns("filt_prot")),
             ### Create an actionButton at the end of filter_protein step to unlock next step
             fluidRow(
               column(12, offset = 0, style = "text-align: right;",
                      h4("Step 3: Move on the next step:"),
                      actionButton(ns("filter_protein_done"), "GO to the Normalization step",
                                   class = "btn-lg btn-success")
               )
             )),

    tabPanel("Normalization",
             value = "glob_normalisation",
             mod_glob_norm_ui(ns("norm")),
             ### Create an actionButton at the end of filter_protein step to unlock next step
             fluidRow(
               column(12, offset = 0, style = "text-align: right;",
                      h4("Step 5: Move on the next step:"),
                      actionButton(ns("normalization_done"), "GO to the Imputation step",
                                   class = "btn-lg btn-success")
               )
             )),

    tabPanel("Imputation",
             value = "glob_imputation",
             mod_glob_imput_ui(ns("imput")),
             ### Create an actionButton at the end of normalisation step to unlock next step
             fluidRow(
               column(12, offset = 0, style = "text-align: right;",
                      h4("Step 3: Move on the next step:"),
                      actionButton(ns("imputation_done"), "GO to the Statistics step",
                                   class = "btn-lg btn-success")
               )
             )),

    tabPanel("Statistics",
             value = "glob_statistics",
             mod_glob_stats_ui(ns("stats")),
             ### Create an actionButton at the end of normalisation step to unlock next step
             fluidRow(
               column(12, offset = 0, style = "text-align: right;",
                      h4("Step 6: Move on to the PHOSPHO analysis:"),
                      actionButton(ns("glob_done"), "GO to the Phospho analysis",
                                   class = "btn-lg btn-success")
               )
             ))

  )
}

#' glob Server Functions
#'
#' @noRd
mod_glob_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    # Step 1 — Create DAlist
    DAlist_obj <- mod_glob_DAlist_server("DAlist")

    # When DAlist creation step is done -> unlock filter_sample step
    observeEvent(input$DAlist_done, {
      shinyjs::enable(selector = "a[data-value='glob_filter_sample']")
      updateTabsetPanel(session, "glob_tabs", selected = "glob_filter_sample")
    })

    # Step 2 — Filter samples
    filt_sample <- mod_glob_filt_sample_server("filt_sample", DAlist_obj)

    # When filter_sample step is done -> unlock filter_protein step
    observeEvent(input$filter_sample_done, {
      shinyjs::enable(selector = "a[data-value='glob_filter_protein']")
      updateTabsetPanel(session, "glob_tabs", selected = "glob_filter_protein")
    })

    # Step 3 — Filter proteins
    filt_prot <- mod_glob_filt_prot_server("filt_prot", filt_sample)

    # When filter_protein step is done -> unlock normalization step
    observeEvent(input$filter_protein_done, {
      shinyjs::enable(selector = "a[data-value='glob_normalisation']")
      updateTabsetPanel(session, "glob_tabs", selected = "glob_normalisation")
    })

    # Step 4 — Normalisation
    norm_data <- mod_glob_norm_server("norm", filt_prot)

    # When Normalisation step is done -> unlock imputation step
    observeEvent(input$normalization_done, {
      shinyjs::enable(selector = "a[data-value='glob_imputation']")
      updateTabsetPanel(session, "glob_tabs", selected = "glob_imputation")
    })

    # Step 5 — Imputation
    imput_data <- mod_glob_imput_server("imput", norm_data)

    # When Imputation step is done -> unlock statistics step
    observeEvent(input$imputation_done, {
      shinyjs::enable(selector = "a[data-value='glob_statistics']")
      updateTabsetPanel(session, "glob_tabs", selected = "glob_statistics")
    })

    # Step 6 — Statistics
    stats_data <- mod_glob_stats_server("stats", imput_data)

    # When ALL GLOBAL STEPS are done -> unlock PHOSPHO analysis
    observeEvent(input$glob_done, {   # glob_done triggered by your button
      session$sendCustomMessage("glob_finished", TRUE)
    })

    # Final reactive result synthesizing all steps
    global_wrap <- reactive({
      req(DAlist_obj(), filt_sample(), filt_prot(), norm_data(), imput_data(), stats_data())

      stats <- stats_data()

      stats$annotation <- harmonize_annotation(
        stats$annotation,
        type = "glob"
      )

      list(
        DAlist = DAlist_obj(),
        filt_sample = filt_sample(),
        filt_prot = filt_prot(),
        norm_data = norm_data(),
        imput_data = imput_data(),
        stats_data = stats
      )
    })

    return(list(
      result = global_wrap
    ))

  })
}

## To be copied in the UI
# mod_glob_ui("glob_1")

## To be copied in the server
# mod_glob_server("glob_1")

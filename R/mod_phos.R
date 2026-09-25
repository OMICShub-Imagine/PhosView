#' phosUI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList

mod_phos_ui <- function(id) {
  ns <- NS(id)

  tabsetPanel(
    id = ns("phos_tabs"),
    tabPanel("1st class phosphosite",
             value = "phos_first_class",
             mod_phos_first_class_ui(ns("first_class")),
             ### Create an actionButton at the end of DAlist creation step to unlock next step
             fluidRow(
               column(12, offset = 0, style = "text-align: right;",
                      h4("Step 5: Move on the next step:"),
                      actionButton(ns("first_class_done"), "GO to the DAlist creation step",
                                   class = "btn-lg btn-success")
               )
             )),

    tabPanel("DAList creation",
             value = "phos_dalist_creation",
             mod_phos_DAlist_ui(
               ns("DAlist")),
             fluidRow(
               column(12, style = "text-align: right;",
                      h4("Step 5: Move on the next step:"),
                      actionButton(ns("DAlist_done"),
                                   "GO to the Filter sample(s) step",
                                   class = "btn-lg btn-success")
               )
             )),

    tabPanel("Filter sample(s)",
             value = "phos_filter_sample",
             mod_phos_filt_sample_ui(ns("filt_sample")),
             ### Create an actionButton at the end of filter_sample(s) step to unlock next step
             fluidRow(
               column(12, offset = 0, style = "text-align: right;",
                      h4("Step 4: Move on the next step:"),
                      actionButton(ns("filter_sample_done"), "GO to the the Filter phosphosite(s) step",
                                   class = "btn-lg btn-success")
               )
             )),

    tabPanel("Filter phosphosite(s)",
             value = "phos_filter_phos",
             mod_phos_filt_phos_ui(ns("filt_phos")),
             ### Create an actionButton at the end of filter_protein step to unlock next step
             fluidRow(
               column(12, offset = 0, style = "text-align: right;",
                      h4("Step 3: Move on the next step:"),
                      actionButton(ns("filter_phos_done"), "GO to the Normalization step",
                                   class = "btn-lg btn-success")
               )
             )),

    tabPanel("Normalisation",
             value = "phos_normalisation",
             mod_phos_norm_ui(ns("norm")),
             ### Create an actionButton at the end of filter_protein step to unlock next step
             fluidRow(
               column(12, offset = 0, style = "text-align: right;",
                      h4("Step 5: Move on the next step:"),
                      actionButton(ns("normalization_done"), "GO to Imputation step",
                                   class = "btn-lg btn-success")
               )
             )),

    tabPanel("Imputation",
             value = "phos_imputation",
             mod_phos_imput_ui(ns("imput")),
             ### Create an actionButton at the end of normalisation step to unlock next step
             fluidRow(
               column(12, offset = 0, style = "text-align: right;",
                      h4("Step 3: Move on the next step:"),
                      actionButton(ns("imputation_done"), "GO to Statistics step",
                                   class = "btn-lg btn-success")
               )
             )),

    tabPanel("Statistics",
             value = "phos_statistics",
             mod_phos_stats_ui(ns("stats")),
             ### Create an actionButton at the end of normalisation step to unlock next step
             fluidRow(
               column(12, offset = 0, style = "text-align: right;",
                      h4("Step 6: Move on to the INTEGRATION analysis:"),
                      actionButton(ns("phos_done"), "GO to INTEGRATION analysis",
                                   class = "btn-lg btn-success")
               )
             ))
  )
}

#' phos Server Functions
#'
#' @noRd
mod_phos_server <- function(id, glob){
  moduleServer(id, function(input, output, session){

    ns <- session$ns

    # Condition
    use_first_class <- reactive({
      req(glob())
      "PG.ProteinGroups" %in% colnames(glob()$stats_data$annotation)
    })

    observe({
      req(use_first_class())  # wait until glob() is available
      if (use_first_class()) {
        shinyjs::hide(selector = paste0("#", NS("DAlist")("browse_wrapper")))
      } else {
        shinyjs::show(selector = paste0("#", NS("DAlist")("browse_wrapper")))
      }
    })



    # Step 1 (optional)
    first_class_obj <- reactiveVal(NULL)

    observe({
      if (use_first_class()) {
        first_class_obj(
          mod_phos_first_class_server("first_class")
        )
      } else {
        # Skip Step 1
        shinyjs::hide(selector = "a[data-value='phos_first_class']")
        shinyjs::enable(selector = "a[data-value='phos_dalist_creation']")

        updateTabsetPanel(session, "phos_tabs",
                          selected = "phos_dalist_creation")
      }
    })

    observeEvent(input$first_class_done, {
      req(use_first_class())
      shinyjs::enable(selector = "a[data-value='phos_dalist_creation']")
      updateTabsetPanel(session, "phos_tabs",
                        selected = "phos_dalist_creation")
    })

    # Step 2 — DAlist
    DAlist_obj <- mod_phos_DAlist_server(
      "DAlist",
      quant_data = reactive({
        if (use_first_class()) {
          req(first_class_obj())
          first_class_obj()$quant_data()
        } else {
          NULL  # will be loaded from browse inside DAlist
        }
      }),
      meta_data = reactive({
        if (use_first_class()) {
          req(first_class_obj())
          first_class_obj()$meta_data()
        } else {
          NULL  # will be loaded from browse inside DAlist
        }
      }),
      use_first_class = use_first_class
    )

    observeEvent(input$DAlist_done, {
      shinyjs::enable(selector = "a[data-value='phos_filter_sample']")
      updateTabsetPanel(session, "phos_tabs",
                        selected = "phos_filter_sample")
    })

    filt_sample <- mod_phos_filt_sample_server("filt_sample", DAlist_obj)

    observeEvent(input$filter_sample_done, {
      shinyjs::enable(selector = "a[data-value='phos_filter_phos']")
      updateTabsetPanel(session, "phos_tabs",
                        selected = "phos_filter_phos")
    })

    filt_phos <- mod_phos_filt_phos_server("filt_phos", filt_sample)

    observeEvent(input$filter_phos_done, {
      shinyjs::enable(selector = "a[data-value='phos_normalisation']")
      updateTabsetPanel(session, "phos_tabs",
                        selected = "phos_normalisation")
    })

    norm_data <- mod_phos_norm_server("norm", filt_phos)

    observeEvent(input$normalization_done, {
      shinyjs::enable(selector = "a[data-value='phos_imputation']")
      updateTabsetPanel(session, "phos_tabs",
                        selected = "phos_imputation")
    })

    imput_data <- mod_phos_imput_server("imput", norm_data)

    observeEvent(input$imputation_done, {
      shinyjs::enable(selector = "a[data-value='phos_statistics']")
      updateTabsetPanel(session, "phos_tabs",
                        selected = "phos_statistics")
    })

    stats_data <- mod_phos_stats_server("stats", imput_data)

    observeEvent(input$phos_done, {
      session$sendCustomMessage("phos_finished", TRUE)
    })

    phospho_wrap <- reactive({
      req(DAlist_obj(), filt_sample(), filt_phos(),
          norm_data(), imput_data(), stats_data())

      stats <- stats_data()

      stats$annotation <- harmonize_annotation(
        stats$annotation,
        type = "phos"
      )

      list(
        DAlist = DAlist_obj(),
        filt_sample = filt_sample(),
        filt_phos = filt_phos(),
        norm_data = norm_data(),
        imput_data = imput_data(),
        stats_data = stats
      )
    })

    return(list(result = phospho_wrap))
  })
}

## To be copied in the UI
# mod_phospho_ui("phospho_1")

## To be copied in the server
# mod_phospho_server("phospho_1")

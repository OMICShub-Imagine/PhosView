#' phos_first_class_browse UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_first_class_browse_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(6,
             wellPanel(
               h4("Step1: Load your quantitative data (TSV)"),
               fileInput(ns("quant_file"), label = NULL, accept = c('.tsv', 'text/tab-separated-values')),
               verbatimTextOutput(ns("upload_status_quant"))
             )
      ),
      column(6,
             wellPanel(
               h4("Step2: Load your experimental metadata (TXT)"),
               fileInput(ns("meta_file"), label = NULL, accept = c('.txt', 'text/plain')),
               verbatimTextOutput(ns("upload_status_meta"))
             )
      )
    )

  )
}

#' phos_first_class_browse Server Functions
#'
#' @noRd
mod_phos_first_class_browse_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns


    # Reactive values to hold uploaded data
    quant_data <- reactiveVal(NULL)
    meta_data <- reactiveVal(NULL)

    # Observe quantitative file upload
    observeEvent(input$quant_file, {
      req(input$quant_file)
      tryCatch({
        quant <- mod_phos_first_class_browse_fct_load(input$quant_file$datapath, data_type = "quant")
        quant_data(quant)
        output$upload_status_quant <- renderText("Quantitative data successfully loaded!")
      }, error = function(e) {
        output$upload_status_quant <- renderText(paste("Error loading quantitative data:", e$message))
      })
    })


    # Observe metadata file upload
    observeEvent(input$meta_file, {
      req(input$meta_file)
      tryCatch({
        meta <- mod_phos_first_class_browse_fct_load(input$meta_file$datapath, data_type = "meta")
        meta_data(meta)
        output$upload_status_meta <- renderText("Metadata successfully loaded!")
      }, error = function(e) {
        output$upload_status_meta <- renderText(paste("Error loading metadata:", e$message))
      })
    })


    # Return reactive values
    return(list(
      quant_data = quant_data,
      meta_data = meta_data
    ))

    return(print(quant_data, meta_date))

  })
}

## To be copied in the UI
# mod_phos_first_class_browse_ui("phos_first_class_browse_1")

## To be copied in the server
# mod_phos_first_class_browse_server("phos_first_class_browse_1")

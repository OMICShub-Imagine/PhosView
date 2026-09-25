#' glob_DAlist_browse UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_glob_DAlist_browse_ui <- function(id) {
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

#' glob_DAlist_browse Server Functions
#'
#' @noRd
mod_glob_DAlist_browse_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    # Reactive values to hold uploaded data
    quant_data <- reactiveVal(NULL)
    meta_data <- reactiveVal(NULL)

    # Observe quantitative file upload
    observeEvent(input$quant_file, {
      req(input$quant_file)
      tryCatch({
        quant <- mod_glob_DAlist_browser_fct_load(input$quant_file$datapath, data_type = "quant")
        quant_data(quant)
        output$upload_status_quant <- renderText("Quantitative data successfully loaded !")
      }, error = function(e) {
        output$upload_status_quant <- renderText(paste("Error loading quantitative data:", e$message))
      })
    })


    # Observe metadata file upload
    observeEvent(input$meta_file, {
      req(input$meta_file)
      tryCatch({
        meta <- mod_glob_DAlist_browser_fct_load(input$meta_file$datapath, data_type = "meta")
        meta_data(meta)
        output$upload_status_meta <- renderText("Metadata successfully loaded !")
      }, error = function(e) {
        output$upload_status_meta <- renderText(paste("Error loading metadata:", e$message))
      })
    })


    # Return reactive values
    return(list(
      quant_data = quant_data,
      meta_data = meta_data
    ))
  })

}

## To be copied in the UI
# mod_glob_DAlist_browse_ui("glob_DAlist_browse_1")

## To be copied in the server
# mod_glob_DAlist_browse_server("glob_DAlist_browse_1")

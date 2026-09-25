#' phos_DAlist_transform UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_DAlist_transform_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             h4("Step 1: Transform your quantitative data and your experimental design file into an integrated R object"),
             actionButton(ns("transform"), "TRANSFORM", class = "btn-lg btn-success"),
             verbatimTextOutput(ns("upload_status"))
      )
    )
  )
}

#' phos_DAlist_transform Server Functions
#'
#' @noRd
mod_phos_DAlist_transform_server <- function(id, quant_data, meta_data){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    # Reactive values to hold uploaded data
    DAlist_obj <- reactiveVal(NULL)

    # Observe quantitative file upload
    observeEvent(input$transform, {
      req(quant_data())
      req(meta_data())
      tryCatch({
        DAlist <- mod_phos_DAlist_transform_fct_DAlist(
          quant_data = quant_data(),
          meta_data = meta_data())
        DAlist_obj(DAlist)
        output$upload_status <- renderText("DAlist object successfully created!")
      }, error = function(e) {
        output$upload_status <- renderText(paste("Error during creation of DAlist object:", e$message))
      })
    })

    # Return the created DAlist so other modules can use it
    return(DAlist_obj)

  })
}

## To be copied in the UI
# mod_phos_DAlist_transform_ui("phos_DAlist_transform_1")

## To be copied in the server
# mod_phos_DAlist_transform_server("phos_DAlist_transform_1")

#' phos_imput_transform UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_imput_transform_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             h4("Step 1: Usually, proteomics workflow implies imputation of missing values.
                Therefore, you can click on the IMPUTATION button just below to impute your missing values"),
             actionButton(ns("imput"), "IMPUTATION", class = "btn-lg btn-success"),
             verbatimTextOutput(ns("upload_status"))
      )
    )
  )
}

#' phos_imput_transform Server Functions
#'
#' @noRd
mod_phos_imput_transform_server <- function(id, DAlist_obj){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    # Reactive values to hold uploaded data
    imput_obj <- reactiveVal(NULL)

    # Observe quantitative file upload
    observeEvent(input$imput, {
      req(DAlist_obj())
      tryCatch({
        imput_data <- mod_phos_imput_transform_fct_imput(
          DAlist_obj = DAlist_obj())
        imput_obj(imput_data)
        output$upload_status <- renderText("DAlist object successfully imputed!")
      }, error = function(e) {
        output$upload_status <- renderText(paste("Error during imputation", e$message))
      })
    })

    # Return the created DAlist so other modules can use it
    return(imput_obj)

  })
}

## To be copied in the UI
# mod_phos_imput_transform_ui("phos_imput_transform_1")

## To be copied in the server
# mod_phos_imput_transform_server("phos_imput_transform_1")

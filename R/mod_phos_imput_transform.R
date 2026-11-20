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
            Therefore, you can click on the IMPUTATION button just below to impute your missing values."),

             p(
               "Missing values are imputed using a Perseus-style down-shifted Gaussian distribution. ",
               "This approach assumes that missing values are predominantly missing not at random (MNAR), ",
               "reflecting low protein abundance below the mass spectrometer detection limit. ",
               "Imputed values are therefore randomly drawn from a Gaussian distribution shifted toward ",
               "lower intensities and narrower than the observed distribution. ",
               "The default parameters (a 1.8 standard-deviation downshift and a distribution width of ",
               "0.3 standard deviations) follow the default parameters in Perseus ",
               "(Tyanova et al., 2016a, 2016b).",
               style = "font-size: 13px; color: #666666; margin-bottom: 15px;"
             ),

             actionButton(
               ns("imput"),
               "IMPUTATION",
               class = "btn-lg btn-success"
             ),

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

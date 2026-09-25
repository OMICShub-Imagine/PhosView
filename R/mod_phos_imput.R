#' phos_imput UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_imput_ui <- function(id) {
  ns <- NS(id)
  tagList(
    mod_phos_imput_transform_ui(ns("imput")),
    h4("Step 2: Generate QC report:"),
    mod_phos_qc_postnorm_ui(ns("qc"))
  )
}

#' phos_imput Server Functions
#'
#' @noRd
mod_phos_imput_server <- function(id, norm_data){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    imput_data <- mod_phos_imput_transform_server(
      id = "imput",
      DAlist_obj = norm_data
    )

    qc <- mod_phos_qc_postnorm_server(
      id = "qc",
      DAlist_obj = imput_data
    )

    return(imput_data)

  })
}

## To be copied in the UI
# mod_phos_imput_ui("phos_imput_1")

## To be copied in the server
# mod_phos_imput_server("phos_imput_1")

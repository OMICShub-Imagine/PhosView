#' glob_imput UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_glob_imput_ui <- function(id) {
  ns <- NS(id)
  tagList(
    mod_glob_imput_transform_ui(ns("imput")),
    h4("Step 2: Generate QC report:"),
    mod_glob_qc_postnorm_ui(ns("qc"))
  )
}

#' glob_imput Server Functions
#'
#' @noRd
mod_glob_imput_server <- function(id, norm_data){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    imput_data <- mod_glob_imput_transform_server(
      id = "imput",
      DAlist_obj = norm_data
    )

    qc <- mod_glob_qc_postnorm_server(
      id = "qc",
      DAlist_obj = imput_data
    )

    return(imput_data)

  })
}

## To be copied in the UI
# mod_glob_imput_ui("glob_imput_1")

## To be copied in the server
# mod_glob_imput_server("glob_imput_1")

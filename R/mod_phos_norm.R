#' phos_norm UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_norm_ui <- function(id) {
  ns <- NS(id)
  tagList(
    mod_phos_norm_report_ui(ns("qc_norm")),
    mod_phos_norm_transform_ui(ns("transform")),
    h4("Step 4: Generate QC report:"),
    mod_phos_qc_postnorm_ui(ns("qc"))

  )
}

#' phos_norm Server Functions
#'
#' @noRd
mod_phos_norm_server <- function(id, filt_phos){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    qc_norm <- mod_phos_norm_report_server(
      id = "qc_norm",
      DAlist_obj = filt_phos
    )

    norm_data <- mod_phos_norm_transform_server(
      id = "transform",
      DAlist_obj = filt_phos
    )

    qc <- mod_phos_qc_postnorm_server(
      id = "qc",
      DAlist_obj = norm_data
    )

    return(norm_data)


  })
}

## To be copied in the UI
# mod_phos_norm_ui("phos_norm_1")

## To be copied in the server
# mod_phos_norm_server("phos_norm_1")

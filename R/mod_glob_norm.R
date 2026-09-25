#' glob_norm UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_glob_norm_ui <- function(id) {
  ns <- NS(id)
  tagList(
    mod_glob_norm_report_ui(ns("qc_norm")),
    mod_glob_norm_transform_ui(ns("transform")),
    h4("Step 4: Generate QC report:"),
    mod_glob_qc_postnorm_ui(ns("qc"))

  )
}

#' glob_norm Server Functions
#'
#' @noRd
mod_glob_norm_server <- function(id, filt_prot){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    qc_norm <- mod_glob_norm_report_server(
      id = "qc_norm",
      DAlist_obj = filt_prot
    )

    norm_data <- mod_glob_norm_transform_server(
      id = "transform",
      DAlist_obj = filt_prot
    )

    qc <- mod_glob_qc_postnorm_server(
      id = "qc",
      DAlist_obj = norm_data
    )

    return(norm_data)


  })
}

## To be copied in the UI
# mod_glob_norm_ui("glob_norm_1")

## To be copied in the server
# mod_glob_norm_server("glob_norm_1")

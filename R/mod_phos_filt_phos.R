#' phos_filt_phos UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_filt_phos_ui <- function(id) {
  ns <- NS(id)
  tagList(
    mod_phos_filt_phos_transform_ui(ns("transform")),
    h4("Step 2: Generate QC report:"),
    mod_phos_qc_prenorm_ui(ns("qc"))
  )
}

#' phos_filt_phos Server Functions
#'
#' @noRd
mod_phos_filt_phos_server <- function(id, filt_sample){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    # Call the transform submodule and pass DAlist_obj
    filt_phos <- mod_phos_filt_phos_transform_server(
      id = "transform",
      DAlist_obj = filt_sample
    )

    # Pass filtered object to QC module if needed
    qc <- mod_phos_qc_prenorm_server(
      id = "qc",
      DAlist_obj = filt_phos
    )

    return(filt_phos)

  })
}

## To be copied in the UI
# mod_phos_filt_phos_ui("phos_filt_phos_1")

## To be copied in the server
# mod_phos_filt_phos_server("phos_filt_phos_1")

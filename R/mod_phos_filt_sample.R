#' phos_filt_sample UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_filt_sample_ui <- function(id) {
  ns <- NS(id)
  tagList(
    mod_phos_filt_sample_transform_ui(ns("transform")),
    h4("Step 3: Generate QC report:"),
    mod_phos_qc_prenorm_ui(ns("qc"))

  )
}

#' phos_filt_sample Server Functions
#'
#' @noRd
mod_phos_filt_sample_server <- function(id, DAlist_obj){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    # Call the transform submodule and pass DAlist_obj
    filt_sample <- mod_phos_filt_sample_transform_server(
      id = "transform",
      DAlist_obj = DAlist_obj
    )

    # Pass filtered object to QC module if needed
    qc <- mod_phos_qc_prenorm_server(
      id = "qc",
      DAlist_obj = filt_sample
    )

    return(filt_sample)

  })
}

## To be copied in the UI
# mod_phos_filt_sample_ui("phos_filt_sample_1")

## To be copied in the server
# mod_phos_filt_sample_server("phos_filt_sample_1")

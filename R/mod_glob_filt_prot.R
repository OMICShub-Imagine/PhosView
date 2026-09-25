#' glob_filt_prot UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_glob_filt_prot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    mod_glob_filt_prot_transform_ui(ns("transform")),
    h4("Step 2: Generate QC report:"),
    mod_glob_qc_prenorm_ui(ns("qc"))
  )
}

#' glob_filt_prot Server Functions
#'
#' @noRd
mod_glob_filt_prot_server <- function(id, filt_sample){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    # Call the transform submodule and pass DAlist_obj
    filt_prot <- mod_glob_filt_prot_transform_server(
      id = "transform",
      DAlist_obj = filt_sample
    )

    # Pass filtered object to QC module if needed
    qc <- mod_glob_qc_prenorm_server(
      id = "qc",
      DAlist_obj = filt_prot
    )

    return(filt_prot)

  })
}

## To be copied in the UI
# mod_glob_filt_prot_ui("glob_filt_prot_1")

## To be copied in the server
# mod_glob_filt_prot_server("glob_filt_prot_1")

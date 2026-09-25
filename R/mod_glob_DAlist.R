#' glob_DAlist UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_glob_DAlist_ui <- function(id) {
  ns <- NS(id)
  tagList(
    mod_glob_DAlist_browse_ui(ns("browse")),
    mod_glob_DAlist_transform_ui(ns("transform")),
    h4("Step 4: Generate QC report:"),
    mod_glob_qc_prenorm_ui(ns("qc"))
  )
}

#' glob_DAlist Server Functions
#'
#' @noRd
mod_glob_DAlist_server <- function(id){
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ### STEP1&2: LOAD QUANT DATA ET ASSOCIATED METADATA
    browse <- mod_glob_DAlist_browse_server("browse")

    ### STEP3: TRANSFORM LOADED DATA INTO AN INTEGRATED R OBJECT
    transform <- mod_glob_DAlist_transform_server(
      id = "transform",
      quant_data = browse$quant_data,
      meta_data  = browse$meta_data
    )

    ### STEP4: PERFORM QUALITY CONTROL OF YOUR DATA
    qc <- mod_glob_qc_prenorm_server(
      id = "qc",
      DAlist_obj = transform
    )

    ### RETURN OBJECT FOR NEXT STEP
    return(transform)

  })
}


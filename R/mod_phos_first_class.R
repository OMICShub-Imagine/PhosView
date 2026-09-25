#' phos_first_class UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_first_class_ui <- function(id) {
  ns <- NS(id)
  tagList(
    mod_phos_first_class_browse_ui(ns("browse")),
    mod_phos_first_class_qc_ui(ns("first_class_qc")),
    mod_phos_first_class_filt_ui(ns("first_class_filt"))
  )
}

#' phos_first_class Server Functions
#'
#' @noRd
mod_phos_first_class_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    browse <- mod_phos_first_class_browse_server("browse")

    first_class_qc <- mod_phos_first_class_qc_server(id = "first_class_qc",
                                                     quant_data = browse$quant_data,
                                                     meta_data = browse$meta_data)

    first_class_phos <- mod_phos_first_class_filt_server(id = "first_class_filt",
                                                     quant_data = browse$quant_data,
                                                     meta_data = browse$meta_data)

    return(
      list(
        quant_data = first_class_phos,
        meta_data  = browse$meta_data
      )
    )

  })
}

## To be copied in the UI
# mod_phos_first_class_ui("phos_first_class_1")

## To be copied in the server
# mod_phos_first_class_server("phos_first_class_1")

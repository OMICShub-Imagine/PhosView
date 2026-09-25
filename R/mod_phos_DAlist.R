#' phos_DAlist UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_DAlist_ui <- function(id, use_first_class) {
  ns <- NS(id)

  tagList(
    # Only include the browse panel if Step 1 is skipped
    uiOutput(ns("browse_ui")),

    # Transform and QC are always visible
    mod_phos_DAlist_transform_ui(ns("transform")),
    h4("Step 2: Generate QC report:"),
    mod_phos_qc_prenorm_ui(ns("qc"))
  )
}

#' phos_DAlist Server Functions
#'
#' @noRd
mod_phos_DAlist_server <- function(id, quant_data = NULL, meta_data = NULL, use_first_class) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Conditionally render the browse UI only if first_class is FALSE
    output$browse_ui <- renderUI({
      if (!use_first_class()) {
        mod_phos_DAlist_browse_ui(ns("browse"))
      } else {
        NULL
      }
    })

    # Load data from browse only if first_class is FALSE
    browse <- mod_phos_DAlist_browse_server("browse")

    # Transform step
    transform <- mod_phos_DAlist_transform_server(
      id = "transform",
      quant_data = reactive({
        if (use_first_class()) {
          req(quant_data())
          quant_data()
        } else {
          req(browse$quant_data())
          browse$quant_data()
        }
      }),
      meta_data = reactive({
        if (use_first_class()) {
          req(meta_data())
          meta_data()
        } else {
          req(browse$meta_data())
          browse$meta_data()
        }
      })
    )

    # QC
    qc <- mod_phos_qc_prenorm_server(
      id = "qc",
      DAlist_obj = transform
    )

    return(transform)
  })
}

## To be copied in the UI
# mod_phos_DAlist_ui("phos_DAlist_1")

## To be copied in the server
# mod_phos_DAlist_server("phos_DAlist_1")

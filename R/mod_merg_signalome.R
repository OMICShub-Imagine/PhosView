#' merg_signalome UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_merg_signalome_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      h4("Kinase substrates heatmap from Up-regulated phosphosites"),
      column(
        width = 12,
        plotly::plotlyOutput(ns("UP_KinSubScore_Heatmap"), height = "1200px")
      ),

      h4("Kinase substrates heatmap from Down-regulated phosphosites"),
      column(
        width = 12,
        plotly::plotlyOutput(ns("DOWN_KinSubScore_Heatmap"), height = "1200px")
      )
    )
  )
}


#' merg_signalome Server Functions
#'
#' @noRd
mod_merg_signalome_server <- function(id, kinsub_results, comparison){
  moduleServer(id, function(input, output, session){

    # UP-regulated heatmap
    output$UP_KinSubScore_Heatmap <- plotly::renderPlotly({
      req(kinsub_results(), comparison())
      renderHeatmap_UP(kinsub_results(), comparison())
    })

    # DOWN-regulated heatmap
    output$DOWN_KinSubScore_Heatmap <- plotly::renderPlotly({
      req(kinsub_results(), comparison())
      renderHeatmap_DOWN(kinsub_results(), comparison())
    })

  })
}


## To be copied in the UI
# mod_merg_signalome_ui("merg_signalome_1")

## To be copied in the server
# mod_merg_signalome_server("merg_signalome_1")

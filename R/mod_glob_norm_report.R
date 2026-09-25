#' glob_norm_report UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_glob_norm_report_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             h4("Step 1: You have the opportunity to generate a QC normalization report that allows you
                 to choose the normalization method that is best for your dataset so you can click on the
                 Generate button just below"),
             actionButton(ns("generate_qc_norm"), "GENERATE", class = "btn-lg btn-success"),
             verbatimTextOutput(ns("qc_status")),

             # Placeholder for all QC content
             uiOutput(ns("qc_norm_content")))

    )

  )
}

#' glob_norm_report Server Functions
#'
#' @noRd
mod_glob_norm_report_server <- function(id, DAlist_obj){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    observeEvent(input$generate_qc_norm, {
      req(DAlist_obj())

      withProgress(message = "Generating QC plots...", value = 0, {
        tryCatch({
          incProgress(0.5, detail = "Creating normalisation report...")

          qc_plots <- mod_glob_norm_report_fct_norm_report2(DAlist_obj(),
                                                            grouping_column = "group")

          incProgress(0.5, detail = "Rendering plots...")

          output$qc_norm_content <- renderUI({
            tagList(
              h4("Plotting metrics to choose best suitable normalization method"),
              fluidRow(
                column(4, plotOutput(ns("a"))),
                column(4, plotOutput(ns("b"))),
                column(4, plotOutput(ns("c")))
              ),
              fluidRow(
                column(4, plotOutput(ns("d"))),
                column(4, plotOutput(ns("e"))),
                column(4, plotOutput(ns("f")))
              ),
              h4("MD plot for each normalization method"),
              fluidRow(
                column(12, plotOutput(ns("g")))
              )
            )
          })

          # Then render individual plots
          output$a <- renderPlot({ qc_plots$a })
          output$b <- renderPlot({ qc_plots$b })
          output$c <- renderPlot({ qc_plots$c })
          output$d <- renderPlot({ qc_plots$d })
          output$e <- renderPlot({ qc_plots$e })
          output$f <- renderPlot({ qc_plots$f })
          output$g <- renderPlot({ qc_plots$g })

          output$qc_status <- renderText("QC plots generated!")

        }, error = function(e) {
          output$qc_status <- renderText(paste("Error generating QC plots:", e$message))
        })
      })
    })

  })
}

## To be copied in the UI
# mod_glob_norm_report_ui("glob_norm_report_1")

## To be copied in the server
# mod_glob_norm_report_server("glob_norm_report_1")

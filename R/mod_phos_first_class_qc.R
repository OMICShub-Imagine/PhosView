#' phos_first_class_qc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_first_class_qc_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Step3: You have now the possibility to generate a qc report based on the quality of phosphosite you want to keep in your analysis.
       Here, qc plots will be reported to investigate the impact of keeping phosphosite having localisation probability >75% in 50% or 70% or 100%
       of the replicates in at least one group. Thanks to that report you will be able to choose the best way to filter your dataset in order to keep
       first class phosphosites."),
    actionButton(ns("generate_qc"), "GENERATE", class = "btn-lg btn-success"),
    verbatimTextOutput(ns("qc_status")),

    # Placeholder for all QC content
    uiOutput(ns("qc_content"))

  )
}

#' phos_first_class_qc Server Functions
#'
#' @noRd
mod_phos_first_class_qc_server <- function(id, quant_data, meta_data){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    observeEvent(input$generate_qc, {

      withProgress(message = "Generating QC plots...", value = 0, {
        tryCatch({
          incProgress(0.2, detail = "Preparing QC data...")
          qc_data <- mod_phos_first_class_qc_utils_data_prep(
            quant_data = quant_data(),
            meta_data = meta_data()
          )

          incProgress(0.6, detail = "Creating QC plots...")
          qc_plots <- mod_phos_first_class_qc_fct_fig(qc_data)

          incProgress(1, detail = "Rendering plots...")

          output$qc_content <- renderUI({
            tagList(
              fluidRow(
                column(6, h4("Number of phosphosites (mean +/- sd) in each filter"),
                       plotOutput(ns("id"))),
                column(6, h4("Ranking plot of phosphosite abundance colored with average localisation probability"),
                       plotOutput(ns("rank_plot")))
              ),
              fluidRow(
                column(6, h4("Mean localisation probability against number of missing LFQ intensities per phosphosite"),
                       plotOutput(ns("missing_vs_loc_boxplot"))),
                column(6, h4("Coefficient of Variation of raw LFQ intensity per phosphosite against Mean localisation probability"),
                       plotOutput(ns("loc_vs_cv_plot")))
              ),
              fluidRow(
                column(6, h4("Distribution  of raw LFQ intensity CVs in each filter"),
                       plotOutput(ns("violin_cv_plot"))),
                column(6, h4("Distribution of mean localisation probability in each filter"),
                       plotOutput(ns("locprob_hist_plot")))
              ),
              fluidRow(
                column(6, h4("CV of localisation probability according to experimental group in each filter"),
                       plotOutput(ns("locprob_per_grp"))),
                column(6, h4("CV of localisation probability according to experimental group and multiplicity in each filter"),
                       plotOutput(ns("locprob_per_M")))
              ),
              fluidRow(
                column(6, h4("Mean LFQ intensity against Mean localisation probability in each filter"),
                       plotOutput(ns("locprob_LFQ_plot"))),
                column(6, h4("CV against Meant LFQ intensity"),
                       plotOutput(ns("CV_LFQ_plot")))
              )
            )
          })

          # Then render individual plots
          output$id <- renderPlot(qc_plots$id)
          output$rank_plot <- renderPlot(qc_plots$rank_plot)
          output$missing_vs_loc_boxplot <- renderPlot(qc_plots$missing_vs_loc_boxplot)
          output$loc_vs_cv_plot <- renderPlot(qc_plots$loc_vs_cv_plot)
          output$violin_cv_plot <- renderPlot(qc_plots$violin_cv_plot)
          output$locprob_hist_plot <- renderPlot(qc_plots$locprob_hist_plot)
          output$locprob_per_grp <- renderPlot(qc_plots$locprob_per_grp)
          output$locprob_per_M <- renderPlot(qc_plots$locprob_per_M)
          output$locprob_LFQ_plot <- renderPlot(qc_plots$locprob_LFQ_plot)
          output$CV_LFQ_plot <- renderPlot(qc_plots$CV_LFQ_plot)


          output$qc_status <- renderText("QC plots generated!")

        }, error = function(e) {
          output$qc_status <- renderText(paste("Error generating QC plots:", e$message))
        })
      })
    })

  })
}

## To be copied in the UI
# mod_phos_first_class_qc_ui("phos_first_class_qc_1")

## To be copied in the server
# mod_phos_first_class_qc_server("phos_first_class_qc_1")

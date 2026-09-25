#' phos_qc_prenorm UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_qc_prenorm_ui <- function(id) {
  ns <- NS(id)
  tagList(
    actionButton(ns("generate_qc"), "GENERATE", class = "btn-lg btn-success"),
    downloadButton(ns("download_qc"), "Download QC Report"),
    verbatimTextOutput(ns("qc_status")),
    # Placeholder for all QC content
    uiOutput(ns("qc_content"))
  )
}

#' phos_qc_prenorm Server Functions
#'
#' @noRd
mod_phos_qc_prenorm_server <- function(id, DAlist_obj){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    qc_report <- reactiveVal(NULL)

    observeEvent(input$generate_qc, {
      req(DAlist_obj())

      withProgress(message = "Generating QC plots...", value = 0, {
        tryCatch({
          incProgress(0.2, detail = "Preparing QC data...")
          qc_data <- mod_phos_qc_prenorm_utils_data_prep(
            quant_data = DAlist_obj()$data,
            meta_data = DAlist_obj()$metadata,
            prot_annot = DAlist_obj()$annotation
          )

          incProgress(0.6, detail = "Creating QC plots...")
          qc_plots <- mod_phos_qc_prenorm_fct_fig(qc_data)

          qc_report(list(
            data = qc_data,
            plots = qc_plots
          ))

          incProgress(1, detail = "Rendering plots...")

          output$qc_content <- renderUI({
            tagList(
              h4("Sample-wise quality control plots"),
              h5("Intensity distribution per sample"),
              fluidRow(
                column(6, plotly::plotlyOutput(ns("s_distrib_int"))),
                column(6, plotly::plotlyOutput(ns("s_bp_int")))
              ),
              h5("Intensity distribution per condition"),
              fluidRow(
                column(6, plotly::plotlyOutput(ns("g_distrib_int"))),
                column(6, plotly::plotlyOutput(ns("g_bp_int")))
              ),
              h5("Number of first-class phosphosites identified per sample and condition"),
              fluidRow(
                column(6, plotly::plotlyOutput(ns("s_id"))),
                column(6, plotly::plotlyOutput(ns("g_id")))
              ),
              h5("Total input intensity per sample and against number of IDs"),
              fluidRow(
                column(6, plotly::plotlyOutput(ns("s_total_int"))),
                column(6, plotly::plotlyOutput(ns("id_vs_int")))
              ),
              h5("Sample correlation and Principal Component Analysis"),
              fluidRow(
                column(6, plotly::plotlyOutput(ns("corr_plot"))),
                column(6, plotOutput(ns("pca_plot")))
              ),
              h5("Heatmap fo the top100 (for each PC1 and PC2) most variable first-class phosphosites"),
              fluidRow(
                column(12, plotly::plotlyOutput(ns("heatmap_pcaload"), width = "90vw", height = "180vh"))
              ),
              h4("Phosphosite-wise quality control plots"),
              h5("Coefficient of Variation against Mean Intensity by condition"),
              fluidRow(
                column(6, plotly::plotlyOutput(ns("g_cv_mean_distrib"))),
                column(6, plotly::plotlyOutput(ns("g_cv_distrib")))
              ),
              h5("Ranking plot"),
              fluidRow(
                column(12, plotly::plotlyOutput(ns("rank_plot")))
              ),
              h5("Valid values distribution across samples and conditions"),
              fluidRow(
                column(6, plotly::plotlyOutput(ns("s_vv_distrib"))),
                column(6, plotly::plotlyOutput(ns("g_vv_distrib")))
              ),
              h5("Heatmap of proteins that are completely absent in one of the conditions"),
              fluidRow(
                column(12, plotly::plotlyOutput(ns("heatmap_on_off"), width = "90vw", height = "90vh"))
              )
            )
          })

          # Then render individual plots
          output$s_distrib_int <- plotly::renderPlotly({ qc_plots$s_distrib_int })
          output$s_bp_int <- plotly::renderPlotly({ qc_plots$s_bp_int })
          output$g_distrib_int <- plotly::renderPlotly({ qc_plots$g_distrib_int })
          output$g_bp_int <- plotly::renderPlotly({ qc_plots$g_bp_int })
          output$s_id <- plotly::renderPlotly({ qc_plots$s_id })
          output$g_id <- plotly::renderPlotly({ qc_plots$g_id })
          output$s_total_int <- plotly::renderPlotly({ qc_plots$s_total_int })
          output$id_vs_int <- plotly::renderPlotly({ qc_plots$id_vs_int })
          output$corr_plot <- plotly::renderPlotly({ qc_plots$corr_plot })
          output$pca_plot <- renderPlot({ qc_plots$pca_plot })
          output$heatmap_pcaload <- plotly::renderPlotly({ qc_plots$heatmap_pcaload })
          output$g_cv_mean_distrib <- plotly::renderPlotly({ qc_plots$g_cv_mean_distrib })
          output$g_cv_distrib <- plotly::renderPlotly({ qc_plots$g_cv_distrib })
          output$rank_plot <- plotly::renderPlotly({ qc_plots$rank_plot })
          output$s_vv_distrib <- plotly::renderPlotly({ qc_plots$s_vv_distrib })
          output$g_vv_distrib <- plotly::renderPlotly({ qc_plots$g_vv_distrib })
          output$heatmap_on_off <- plotly::renderPlotly({ qc_plots$heatmap_on_off })

          output$qc_status <- renderText("QC plots generated!")

        }, error = function(e) {
          output$qc_status <- renderText(paste("Error generating QC plots:", e$message))
        })
      })
    })

    output$download_qc <- downloadHandler(
      filename = function() {
        paste0("qc_report_",
               format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
               ".html")
      },
      content = function(file) {

        req(qc_report())

        temp_dir <- tempdir()

        # Save QC object
        qc_path <- file.path(temp_dir, "qc.rds")
        saveRDS(qc_report(), qc_path)

        # Copy template
        temp_qmd <- file.path(temp_dir, "QC_report_pre_norm.qmd")
        template <- system.file("quarto/QC_report_pre_norm.qmd", package = "PhosView")
        file.copy(template, temp_qmd, overwrite = TRUE)

        # Define output filename (NO PATH!)
        out_file <- "qc_report.html"

        # Render in tempdir
        quarto::quarto_render(
          input = temp_qmd,
          output_file = out_file,
          execute_params = list(
            input_file = qc_path
          )
        )

        # Copy result to Shiny download path
        file.copy(file.path(temp_dir, out_file), file, overwrite = TRUE)
      }
    )

  })
}

## To be copied in the UI
# mod_phos_qc_prenorm_ui("phos_qc_prenorm_1")

## To be copied in the server
# mod_phos_qc_prenorm_server("phos_qc_prenorm_1")

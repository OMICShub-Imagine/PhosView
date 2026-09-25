#' export UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_export_ui <- function(id) {
  ns <- NS(id)

  tagList(

    tags$h3("Statistical analyses export:"),

    ## ---- GLOBAL ANALYSIS EXPORT ----
    tags$h4(
      "Export all quantitative tables from each step of the global analysis"
    ),

    downloadButton(
      outputId = ns("download_global_xls"),
      label = tagList(icon("file-excel"), "Download Global analysis results (.xls)")
    ),

    downloadButton(
      outputId = ns("download_global_rds"),
      label = tagList(icon("database"), "Download Global analysis object (.rds)")
    ),

    tags$hr(),

    ## ---- PHOSPHO ANALYSIS EXPORT ----
    tags$h4(
      "Export all quantitative tables from each step of the phospho analysis"
    ),

    downloadButton(
      outputId = ns("download_phospho_xls"),
      label = tagList(icon("file-excel"), "Download Phospho analysis results (.xls)")
    ),

    downloadButton(
      outputId = ns("download_phospho_rds"),
      label = tagList(icon("database"), "Download Phospho analysis object (.rds)")
    ),

    ## ---- FUNCTIONAL ANALYSIS EXPORT ----

    tags$h3("Functional analyses export:"),

    ### ORA RESULTS EXPORT BUTTONS (PHOS UP/DOWN and GLOB UP/DOWN)
    tags$h4("ORA results:"),
    downloadButton(
      ns("download_ora_xls"),
      label = tagList(icon("file-excel"),
                      "Download ORA analysis results (.xls)")
      ),
    downloadButton(
      ns("download_ora_rds"),
      label = tagList(icon("database"),
                      "Download ORA analysis object (.rds)")
      ),

    ### KINASE-SUBSTRATE RESULTS EXPORT BUTTONS (PHOS UP/DOWN and GLOB UP/DOWN)
    tags$h4("Kinase-substrate results:"),
    downloadButton(
      ns("download_kinsub_xls"),
      label = tagList(icon("file-excel"),
                      "Download kinase/substrate results (.xls)")
    ),

    downloadButton(
      ns("download_kinsub_rds"),
      label = tagList(icon("database"),
                      "Download kinase/substrate results (.rds)")
    ),


    ### PATHVIEW/KEGG RESULTS EXPORT BUTTONS (PHOS UP/DOWN and GLOB UP/DOWN)
    tags$h4("KEGG-Pathview results:"),
    downloadButton(
      ns("download_pathview_xls"),
      label = tagList(icon("file-excel"),
                      "Download Pathview analysis results (.xls)")
      ),
    downloadButton(
      ns("download_pathview_rds"),
      label = tagList(icon("database"),
                      "Download Pathview analysis results (.rds)")
      ),

    ### COMPLETE SESSION EXPORT

    tags$hr(),

    tags$h3("Session export"),

    tags$p(
      "Save the complete PhosView session to reload it later."
    ),

    downloadButton(
      ns("download_session_rds"),
      label = tagList(
        icon("save"),
        "Download complete session (.rds)"
      )
    )
  )
}

#' export Server Functions
#'
#' @noRd
mod_export_server <- function(id, glob, phos, ora, kinsub, pathview, taxo) {
  moduleServer(id, function(input, output, session) {

    ora_res      <- reactive({ req(ora()) })
    kinsub_res   <- reactive({ req(kinsub()) })
    pathview_res <- reactive({ req(pathview()) })

    ## ===============================
    ## GLOBAL STATS EXPORT
    ## ===============================
    output$download_global_xls <- downloadHandler(

      filename = function() {
        paste0(
          "Global_DA_results_",
          format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
          ".xlsx"
        )
      },

      content = function(file) {

        res <- req(glob$result())
        wb  <- openxlsx::createWorkbook()

        for (nm in names(res)) {

          if (is.null(res[[nm]]$data)) next

          openxlsx::addWorksheet(wb, nm)

          if (nm == "stats_data") {
            df <- combine_stats(res[[nm]])
          } else {
            df <- combine_data_annotation(res[[nm]])
          }

          ## Write table starting at row 2
          openxlsx::writeData(wb, nm, df, startRow = 2)

          ## Column groups
          annot_cols <- get_annot_cols(res[[nm]], phospho = FALSE)
          quant_cols <- get_quant_cols(df, res[[nm]])

          stats_groups <- NULL
          if (nm == "stats_data" && !is.null(res[[nm]]$results)) {
            stats_groups <- get_stats_groups(df)
          }

          ## Add merged header row
          add_grouped_headers(
            wb = wb,
            sheet = nm,
            df = df,
            annot_cols = annot_cols,
            quant_cols = quant_cols,
            stats_groups = stats_groups
          )

          if (nm %in% c("imput_data", "stats_data") && "norm_data" %in% names(res)) {

            norm_mat  <- res$norm_data$data
            imput_mat <- res$imput_data$data
            imputed   <- get_imputed_matrix(norm_mat, imput_mat)

            data_cols <- which(colnames(df) %in% colnames(imput_mat))

            grey_style <- openxlsx::createStyle(fgFill = "#D9D9D9")

            for (j in seq_along(data_cols)) {
              rows_to_style <- which(imputed[, j])
              if (length(rows_to_style) > 0) {
                openxlsx::addStyle(
                  wb,
                  sheet = nm,
                  style = grey_style,
                  rows  = rows_to_style + 2,
                  cols  = data_cols[j],
                  stack = TRUE
                )
              }
            }
          }
        }

        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )

    output$download_global_rds <- downloadHandler(

      filename = function() {
        paste0(
          "Global_DA_object_",
          format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
          ".rds"
        )
      },

      content = function(file) {
        res <- req(glob$result())
        saveRDS(res, file)
      }
    )

    ## ===============================
    ## PHOSPHO STATS EXPORT
    ## ===============================
    output$download_phospho_xls <- downloadHandler(

      filename = function() {
        paste0(
          "Phospho_DA_results_",
          format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
          ".xlsx"
        )
      },

      content = function(file) {

        res <- req(phos$result())
        wb  <- openxlsx::createWorkbook()

        for (nm in names(res)) {

          if (is.null(res[[nm]]$data)) next

          openxlsx::addWorksheet(wb, nm)

          if (nm == "stats_data" && !is.null(res[[nm]]$results)) {
            df <- combine_stats_phos(res[[nm]])
          } else {
            df <- combine_data_annotation_phos(res[[nm]])
          }

          ## Write table starting at row 2
          openxlsx::writeData(wb, nm, df, startRow = 2)

          ## Column groups
          annot_cols <- get_annot_cols(res[[nm]], phospho = TRUE)
          quant_cols <- get_quant_cols(df, res[[nm]])

          stats_groups <- NULL
          if (nm == "stats_data" && !is.null(res[[nm]]$results)) {
            stats_groups <- get_stats_groups(df)
          }

          ## Add merged header row
          add_grouped_headers(
            wb = wb,
            sheet = nm,
            df = df,
            annot_cols = annot_cols,
            quant_cols = quant_cols,
            stats_groups = stats_groups
          )

          if (nm %in% c("imput_data", "stats_data") && "norm_data" %in% names(res)) {

            norm_mat  <- res$norm_data$data
            imput_mat <- res$imput_data$data
            imputed   <- get_imputed_matrix(norm_mat, imput_mat)

            data_cols <- which(colnames(df) %in% colnames(imput_mat))

            grey_style <- openxlsx::createStyle(fgFill = "#D9D9D9")

            for (j in seq_along(data_cols)) {
              rows_to_style <- which(imputed[, j])
              if (length(rows_to_style) > 0) {
                openxlsx::addStyle(
                  wb,
                  sheet = nm,
                  style = grey_style,
                  rows  = rows_to_style + 2,
                  cols  = data_cols[j],
                  stack = TRUE
                )
              }
            }
          }
        }

        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )

    output$download_phospho_rds <- downloadHandler(

      filename = function() {
        paste0(
          "Phospho_DA_object_",
          format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
          ".rds"
        )
      },

      content = function(file) {
        res <- req(phos$result())
        saveRDS(res, file)
      }
    )

    ## ===============================
    ## ORA EXPORT
    ## ===============================

    output$download_ora_xls <- downloadHandler(

      filename = function() {
        paste0(
          "ORA_results_",
          format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
          ".xlsx"
        )
      },

      content = function(file) {

        wb <- openxlsx::createWorkbook()
        ora <- req(ora_res())

        for (mode in c("global", "phospho")) {

          tables <- extract_ora_tables(ora[[mode]])

          for (nm in names(tables)) {
            sheet <- paste(mode, nm, sep = "_")
            openxlsx::addWorksheet(wb, sheet)
            openxlsx::writeData(wb, sheet, tables[[nm]])
          }
        }

        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )

    output$download_ora_rds <- downloadHandler(

      filename = function() {
        paste0(
          "ORA_results_",
          format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
          ".rds"
        )
      },

      content = function(file) {
        saveRDS(req(ora_res()), file)
      }
    )

    ## ===============================
    ## KINSUB EXPORT
    ## ===============================

    output$download_kinsub_xls <- downloadHandler(

      filename = function() {
        paste0(
          "Kinase_Substrate_results_",
          format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
          ".xlsx"
        )
      },

      content = function(file) {

        wb <- openxlsx::createWorkbook()
        sig <- req(kinsub_res())

        tables <- extract_signalome_tables(sig)

        for (nm in names(tables)) {
          openxlsx::addWorksheet(wb, nm)
          openxlsx::writeData(wb, nm, tables[[nm]], rowNames = TRUE)
        }

        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )


    output$download_kinsub_rds <- downloadHandler(

      filename = function() {
        paste0(
          "Kinase_Substrate_results_",
          format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
          ".rds"
        )
      },

      content = function(file) {
        saveRDS(req(kinsub_res()), file)
      }
    )

    ## ===============================
    ## PATHVIEW/KEGG EXPORT
    ## ===============================

    output$download_pathview_xls <- downloadHandler(

      filename = function() {
        paste0(
          "Pathview_KEGG_results_",
          format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
          ".xlsx"
        )
      },

      content = function(file) {

        wb <- openxlsx::createWorkbook()
        pv <- req(pathview_res())

        ## Enrichment tables
        for (cmp in names(pv$enrichment)) {

          kegg_res <- pv$enrichment[[cmp]]$KEGG@result
          sheet <- paste0("KEGG_", cmp)

          openxlsx::addWorksheet(wb, sheet)
          openxlsx::writeData(wb, sheet, kegg_res)
        }

        ## Final input used for pathview
        openxlsx::addWorksheet(wb, "Pathview_input")
        openxlsx::writeData(wb, "Pathview_input", pv$final_input)

        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )

    output$download_pathview_rds <- downloadHandler(

      filename = function() {
        paste0(
          "Pathview_KEGG_results_",
          format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
          ".rds"
        )
      },

      content = function(file) {
        saveRDS(req(pathview_res()), file)
      }
    )

    output$download_session_rds <- downloadHandler(

      filename = function() {
        paste0(
          "PhosView_session_",
          format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
          ".rds"
        )
      },

      content = function(file) {

        snapshot <- list(

          glob = req(glob$result()),

          phos = req(phos$result()),

          ora = req(ora_res()),

          kinsub = req(kinsub_res()),

          pathview = req(pathview_res()),

          taxo = input$taxo

        )

        saveRDS(snapshot, file)

      }
    )

  })
}


## To be copied in the UI
# mod_export_ui("export_1")

## To be copied in the server
# mod_export_server("export_1")

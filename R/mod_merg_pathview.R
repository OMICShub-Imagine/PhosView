#' merg_pathview UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom utils data head read.table str
#'
#' @importFrom shiny NS tagList
mod_merg_pathview_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      h4("KEGG ORA on Dysregulated Kinase/Substrates"),
      column(12, DT::dataTableOutput(ns("Kinase_Substrate_KEGG_Table"))),
      hr(),
      h4("Table of phosphosite (substrate) and Kinase retained"),
      column(12, DT::dataTableOutput(ns("Kinase_Substrate_Table"))),
      hr(),
      h4("Pathviewer of selected pathway"),
      column(12, imageOutput(ns("Kinase_Substrate_KEGG_Pathway")))

      )
    )
}

#' merg_pathview Server Functions
#'
#' @noRd
mod_merg_pathview_server <- function(id,
                                     kinsub_ora_selected_results,
                                     kinsub_selected_results,
                                     comparison,
                                     kegg_selector,
                                     taxo_kegg,
                                     AnnotDb,
                                     pathview_dir
                                     ){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    # ---- KEGG Table ----
    output$Kinase_Substrate_KEGG_Table <- DT::renderDataTable({
      res <- kinsub_ora_selected_results()

      if (is.null(res)) {
        validate(need(FALSE, "No ORA results"))
      }

      if (is.null(res$KEGG)) {
        validate(need(FALSE, "No KEGG enrichment results"))
      }

      if (is.null(res$KEGG@result) || nrow(res$KEGG@result) == 0) {
        validate(need(FALSE, "Empty KEGG result table"))
      }

      DT::datatable(as.data.frame(res$KEGG@result),
                    class = "cell-border compact hover order-column",
                    rownames = FALSE,
                    editable = "row",
                    extensions = c("Scroller", "FixedColumns", "Buttons"),
                    selection = "none",
                    filter = list(position = 'top', clear = FALSE),
                    callback = htmlwidgets::JS(
                      "table.columns.adjust().draw();"
                    ),
                    options = list(
                      columnDefs = list(
                        list(
                          targets = c("Description"),
                          render = htmlwidgets::JS(
                            "function(data, type, row, meta) {",
                            "  if (type === 'display' && data !== null && data !== undefined && data !== '') {",
                            "    return data.length > 25 ? '<span title=\"' + data + '\">' + data.substr(0, 25) + '...</span>' : data;",
                            "  } else {",
                            "    return '';",
                            "  }",
                            "}"
                          )
                        )),
                      autoWidth = TRUE,
                      dom = "frtipB",
                      scrollY = 200, scrollX = 400, scroller = TRUE,
                      fixedColumns = list(leftColumns = 4),
                      select = list(style = 'os', items = 'row'),
                      buttons = c('copy', 'csv', 'excel', 'pdf', 'print',
                                  'colvis','selectAll', 'selectNone',
                                  'selectRows', 'selectColumns', 'selectCells')))
    })

    # ---- Kinase/Substrate retained for each comparison (Table) ----
    output$Kinase_Substrate_Table <- DT::renderDataTable({
      res <- kinsub_selected_results()

      if (is.null(res)) {
        validate(need(FALSE, "No ORA results"))
      }


      DT::datatable(as.data.frame(res),
                    class = "cell-border compact hover order-column",
                    rownames = FALSE,
                    editable = "row",
                    extensions = c("Scroller", "FixedColumns", "Buttons"),
                    selection = "none",
                    filter = list(position = 'top', clear = FALSE),
                    callback = htmlwidgets::JS(
                      "table.columns.adjust().draw();"
                    ),
                    options = list(
                      columnDefs = list(
                        list(
                          targets = c("gene"),
                          render = htmlwidgets::JS(
                            "function(data, type, row, meta) {",
                            "  if (type === 'display' && data !== null && data !== undefined && data !== '') {",
                            "    return data.length > 25 ? '<span title=\"' + data + '\">' + data.substr(0, 25) + '...</span>' : data;",
                            "  } else {",
                            "    return '';",
                            "  }",
                            "}"
                          )
                        )),
                      autoWidth = FALSE,
                      dom = "frtipB",
                      scrollY = 200, scrollX = 400, scroller = TRUE,
                      select = list(style = 'os', items = 'row'),
                      buttons = c('copy', 'csv', 'excel', 'pdf', 'print',
                                  'colvis','selectAll', 'selectNone',
                                  'selectRows', 'selectColumns', 'selectCells'))
                    )
    })


    #==============================================================
    # Pathview rendering from KEGG selector
    #==============================================================

    # ---- Build geneList for Pathview ----
    geneList_pathview <- reactive({

      message("\n=== DEBUG geneList_pathview ===")

      df <- kinsub_selected_results()

      message("KinSub results NULL: ", is.null(df))

      if (!is.null(df)) {
        message("KinSub rows: ", nrow(df))
        message("KinSub columns: ", paste(names(df), collapse = ", "))
      }

      message("Retrieving annotation database...")

      annot <- AnnotDb()

      message("Annotation database retrieved")
      message("Annotation database NULL: ", is.null(annot))
      message("Annotation database class: ", paste(class(annot), collapse = ", "))

      message("Annotation database NULL: ", is.null(annot))

      req(!is.null(df), nrow(df) > 0, !is.null(annot))

      message("KinSub results and annotation database validated")

      print(is.null(kinsub_selected_results()))

      if (!is.null(kinsub_selected_results()))
        print(names(kinsub_selected_results()))
      message("geneList_pathview() entered")

      df <- kinsub_selected_results()

      if (is.null(df) || nrow(df) == 0 || all(is.na(df$gene))) {
        message("Empty kinsub_selected_results() : skipping Pathview")
        return(numeric(0))
      }

      genes <- unique(df$gene)
      genes <- unique(sapply(strsplit(genes, ";"), `[`, 1))
      cat("Unique gene symbols after split:", length(genes), "\n")
      cat("Example genes:", paste(head(genes, 10), collapse = ", "), "\n")

      if (length(genes) == 0) {
        message("No valid genes : skipping Pathview")
        return(numeric(0))
      }

      entrez <- AnnotationDbi::mapIds(
        annot,
        keys      = genes,
        column    = "ENTREZID",
        keytype   = "SYMBOL",
        multiVals = "first"
      )

      cat("Mapped ENTREZ IDs:", sum(!is.na(entrez)), "/", length(entrez), "\n")
      cat("Example ENTREZ mapping:\n")
      print(head(entrez, 10))

      df <- df %>%
        dplyr::mutate(gene = strsplit(gene, ";")) %>%
        tidyr::unnest(gene) %>%
        dplyr::filter(gene %in% names(entrez)) %>%
        dplyr::distinct(gene, .keep_all = TRUE)

      cat("Rows after gene filtering:", nrow(df), "\n")
      cat("Structrure df:", str(df), "\n")

      # ---- logFC diagnostics ----
      cat("Non-NA logFC values:", sum(!is.na(df$logFC)), "/", nrow(df), "\n")
      cat("logFC summary:\n")
      print(summary(df$logFC))


      # ---- Build geneList ----
      mapped <- entrez[df$gene]

      geneList <- df$logFC
      names(geneList) <- mapped

      geneList <- geneList[!is.na(names(geneList))]
      geneList <- geneList[!is.na(geneList)]
      geneList[is.na(geneList)] <- 0

      cat("Final geneList length:", length(geneList), "\n")
      cat("Example geneList entries:\n")
      print(geneList)

      message("geneList built")
      return(geneList)

    })


    # ---- Generate Pathview PNG when KEGG changes ----

    pathview_png <- reactive({

      message("=== REGENERATING KEGG PATHWAY ===")

      # 1. Rebuild geneList from saved KinSub results
      geneList <- geneList_pathview()

      validate(
        need(!is.null(geneList), "No genes available"),
        need(length(geneList) > 0, "No enriched phosphosites")
      )

      # 2. Retrieve selected KEGG pathway
      pid <- kegg_selector()

      validate(
        need(!is.null(pid), "No KEGG pathway selected"),
        need(length(pid) == 1, "Invalid KEGG selection"),
        need(!is.na(pid), "KEGG pathway is NA"),
        need(nchar(trimws(pid)) > 0, "Empty KEGG pathway ID")
      )

      # 3. Retrieve taxonomy and output directory
      species <- taxo_kegg()
      out_dir <- pathview_dir

      req(species, out_dir)

      # Ensure directory exists
      if (!dir.exists(out_dir)) {
        dir.create(
          out_dir,
          recursive = TRUE,
          showWarnings = FALSE
        )
      }

      # 4. Generate a fresh pathway image
      message("Pathway ID: ", pid)
      message("Species: ", species)
      message("Genes: ", length(geneList))
      message("Output directory: ", out_dir)

      img <- generate_pathview_png(
        geneList   = geneList,
        pathway_id = pid,
        taxo_kegg  = species,
        out_dir    = out_dir,
        suffix = paste0(
          comparison(),
          "_pathview_",
          basename(tempfile(pattern = "generation_"))
        )
      )

      # 5. Validate generated image
      validate(
        need(!is.null(img), "Pathview generation failed"),
        need(file.exists(img), "Generated Pathview PNG not found")
      )

      message("Pathway successfully generated: ", img)

      return(img)

    })



    # ---- Render Pathview image ----
    output$Kinase_Substrate_KEGG_Pathway <- renderImage({

      req(pathview_png())

      img <- pathview_png()

      validate(
        need(!is.null(img), "No Pathview image available"),
        need(file.exists(img), "Pathview file missing")
      )

      list(
        src = normalizePath(img, winslash = "/", mustWork = TRUE),
        contentType = "image/png"
      )
    }, deleteFile = FALSE)

  })

}

## To be copied in the UI
# mod_merg_pathview_ui("merg_pathview_1")

## To be copied in the server
# mod_merg_pathview_server("merg_pathview_1")

#' merg_dual_view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_merg_dual_view_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12, DT::dataTableOutput(ns("GlobalTable")))
    ),

    hr(),

    fluidRow(
      column(6, plotly::plotlyOutput(ns("input_VP"))),
      column(6, plotly::plotlyOutput(ns("phospho_VP")))
    ),

    hr(),

    fluidRow(
      column(6, plotly::plotlyOutput(ns("input_PP"))),
      column(6, plotly::plotlyOutput(ns("phospho_PP")))
    ),

    hr(),

    fluidRow(
      column(12, plotly::plotlyOutput(ns("FC_plot")))
    )

  )
}

#' merg_dual_view Server Functions
#'
#' @noRd
mod_merg_dual_view_server <- function(id, glob, phos, comparison, user_gene_list, restored_session = reactive(NULL)){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    glob_data <- reactive({

      sess <- restored_session()

      if (!is.null(sess)) {
        return(sess$glob)
      }

      req(glob$result())
      glob$result()

    })

    phos_data <- reactive({

      sess <- restored_session()

      if (!is.null(sess)) {
        return(sess$phos)
      }

      req(phos$result())
      phos$result()

    })

    # Reactive expression for preprocessed data
    data_final <- reactive({
      req(glob_data(), phos_data())

      if("PG.ProteinGroups" %in% colnames(glob_data()$stats_data$annotation)) {
        # Spectronaut
        mod_merg_dual_view_utils_preprocess(
          glob_stats          = glob_data()$stats_data,
          phos_stats          = phos_data()$stats_data,
          protein_id_col_glob = "PG.ProteinGroups",
          protein_id_col_phos = "PG.ProteinGroups",
          gene_col_glob       = "PG.Genes",
          gene_col_phos       = "PG.Genes",
          name_col_glob       = "PG.ProteinNames",
          name_col_phos       = "PG.ProteinNames",
          desc_col_glob       = "PG.ProteinDescriptions",
          desc_col_phos       = "PG.ProteinDescriptions",
          phospho_site_col    = "uniprot_id"
        )
      } else {
        # DIANN
        mod_merg_dual_view_utils_preprocess(
          glob_stats          = glob_data()$stats_data,
          phos_stats          = phos_data()$stats_data,
          protein_id_col_glob = "Protein.Group",
          protein_id_col_phos = "Protein",
          gene_col_glob       = "Genes",
          gene_col_phos       = "Gene.Names",
          name_col_glob       = "Protein.Names",
          name_col_phos       = "Protein.Names",
          desc_col_glob       = "First.Protein.Description",
          desc_col_phos       =  NULL,
          phospho_site_col    = "uniprot_id"
        )
      }
    })


    output$GlobalTable = DT::renderDataTable({
      req(comparison())
      DT::datatable(data_final()$subset_data_final,
                    class = "cell-border compact hover order-column",
                    rownames = FALSE,
                    editable = "row",
                    extensions = c("Scroller", "FixedColumns", "Buttons"),
                    selection = "none",
                    filter = list(position = 'top', clear = FALSE),
                    options = list(
                      columnDefs = list(
                        list(
                          targets = grep("Description", colnames(data_final()$subset_data_final), ignore.case = TRUE, value = TRUE),
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
                      # scroll :
                      scrollY = 200, scrollX = 400, scroller = TRUE,
                      # fixer les colonnes :
                      fixedColumns = list(leftColumns = 2),
                      # selection :
                      select = list(style = 'os', items = 'row'),
                      buttons = c(
                        # enregistrements
                        'copy', 'csv', 'excel', 'pdf', 'print',
                        # visualisation des colonnes
                        'colvis',
                        # selection des elements
                        'selectAll', 'selectNone', 'selectRows', 'selectColumns', 'selectCells'))) %>%
        DT::formatStyle(c(paste0(comparison(), ".sig.FDR_GlobalAnalysis"),
                      paste0(comparison(), ".sig.FDR_PhosphoAnalysis")),
                    backgroundColor  = DT::styleEqual(c(-1, 0, 1), c("#5555FF", "#AAAAAA", "#FF5555"))
        )
    })

    col_config <- reactive({
      req(glob_data())

      if ("PG.ProteinGroups" %in% colnames(glob_data()$stats_data$annotation)) {
        # Spectronaut
        list(
          gene_col_glob = "PG.Genes",
          gene_col_phos = "PG.Genes",
          protein_id_col_glob = "PG.ProteinGroups",
          protein_id_col_phos = "PG.ProteinGroups",
          desc_col_phos = "PG.ProteinDescriptions",
          phospho_site_col = "uniprot_id"
        )
      } else {
        # DIANN
        list(
          gene_col_glob = "Genes",
          gene_col_phos = "Gene.Names",
          protein_id_col_glob = "Protein.Group",
          protein_id_col_phos = "Protein",
          desc_col_phos = NULL,
          phospho_site_col = "uniprot_id"
        )
      }
    })


    # Volcano plots
    output$input_VP <- plotly::renderPlotly({
      req(comparison())
      renderVolcano_Glob(
        data = data_final()$subset_data_final,
        comparison = comparison(),
        user_gene_list = user_gene_list()
      )
    })

    output$phospho_VP <- plotly::renderPlotly({
      req(comparison())
      renderVolcano_Phos(
        data = data_final()$subset_data_final,
        comparison = comparison(),
        user_gene_list = user_gene_list()
      )
    })


    # Profile plots
    output$input_PP <- plotly::renderPlotly({
      cfg <- col_config()

      renderPP_Glob(
        data = data_final()$global_data_final,
        metadata = glob_data()$stats_data$metadata,
        user_gene_list = user_gene_list(),
        gene_col_glob = cfg$gene_col_glob,
        protein_id_col_glob = cfg$protein_id_col_glob
      )
    })

    output$phospho_PP <- plotly::renderPlotly({
      cfg <- col_config()

      renderPP_Phos(
        data = data_final()$phospho_data_final,
        metadata = phos_data()$stats_data$metadata,
        user_gene_list = user_gene_list(),
        gene_col_phos = cfg$gene_col_phos,
        protein_id_col_phos = cfg$protein_id_col_phos,
        phospho_site_col = cfg$phospho_site_col
      )
    })

    # FC plot
    output$FC_plot <- plotly::renderPlotly({
      req(comparison())
      cfg <- col_config()

      renderFC(
        global_data = data_final()$global_data_final,
        phospho_data = data_final()$phospho_data_final,
        comparison = comparison(),
        user_gene_list = user_gene_list(),
        gene_col_glob = cfg$gene_col_glob,
        gene_col_phos = cfg$gene_col_phos,
        protein_id_col_glob = cfg$protein_id_col_glob,
        protein_id_col_phos = cfg$protein_id_col_phos,
        desc_col_phos = cfg$desc_col_phos,
        phospho_site_col = cfg$phospho_site_col
      )
    })

  })
}

## To be copied in the UI
# mod_merg_dual_view_ui("dual_view_1")

## To be copied in the server
# mod_merg_dual_view_server("dual_view_1")

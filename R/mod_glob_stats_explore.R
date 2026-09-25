#' glob_stats_explore UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_glob_stats_explore_ui <- function(id) {
  ns <- NS(id)

  tagList(
    h4("Step 5: Explore statistical results"),

    fluidRow(
      column(
        4,
        uiOutput(ns("comparison_ui"))
      )
    ),

    fluidRow(
      column(12, DT::dataTableOutput(ns("GlobalTable"))),
      uiOutput(ns("stats_vizu_content"))
    ),

    fluidRow(
      column(4,
             selectizeInput(
               ns("user_gene_list"),
               label = "Select gene/protein to highlight:",
               choices = NULL,
               multiple = FALSE
             ))
    ),

    fluidRow(
      column(6, plotly::plotlyOutput(ns("input_VP"), height = "400px")),
      column(6, plotly::plotlyOutput(ns("input_PP"), height = "400px"))
    )
  )
}

#' glob_stats_explore Server Functions
#'
#' @noRd
mod_glob_stats_explore_server <- function(id, DAlist_obj){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    col_info <- reactive({
      req(DAlist_obj())

      if ("Protein.Group" %in% colnames(DAlist_obj()$annotation)) {
        list(
          col_pg   = "Protein.Group",
          col_gene = "Genes",
          col_desc = "First.Protein.Description"
        )
      } else if ("PG.ProteinGroups" %in% colnames(DAlist_obj()$annotation)) {
        list(
          col_pg   = "PG.ProteinGroups",
          col_gene = "PG.Genes",
          col_desc = "PG.ProteinDescriptions"
        )
      }
    })

    # ----- CREATE comparison selector UI -----
    output$comparison_ui <- renderUI({
      req(DAlist_obj())
      comps <- names(DAlist_obj()$results)

      radioButtons(
        ns("comparison"),
        label = "Select comparison:",
        choices = comps,
        selected = comps[1],
        inline = TRUE
      )
    })

    # ----- UPDATE gene selector -----
    observe({
      req(DAlist_obj(), col_info())

      updateSelectizeInput(
        session,
        "user_gene_list",
        choices = sort(DAlist_obj()$annotation[[col_info()$col_gene]]),
        server = TRUE
      )
    })

    # ------------------------- DATATABLE -----------------------------
    output$GlobalTable <- DT::renderDataTable({
      req(DAlist_obj(), input$comparison)

      df <- DAlist_obj()$results[[input$comparison]] %>% tibble::rownames_to_column("uniprot_id") %>% dplyr::left_join(DAlist_obj()$annotation, by = "uniprot_id")
      req(df)

      DT::datatable(df,
                    class = "cell-border compact hover order-column",
                    rownames = FALSE,
                    editable = "row",
                    extensions = c("Scroller", "FixedColumns", "Buttons"),
                    selection = "none",
                    filter = list(position = 'top', clear = FALSE),
                    options = list(
                      columnDefs = list(
                        list(
                          targets = col_info()$col_desc,
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
        DT::formatStyle("sig.FDR",
                        backgroundColor  = DT::styleEqual(c(-1, 0, 1), c("#5555FF", "#AAAAAA", "#FF5555"))
        )

    })

    # ---------------------- CREATE VP AND PP OBJECT -----------------------------

    stats_plot <- reactive({
      req(DAlist_obj(), input$comparison)

      mod_glob_stats_explore_fct_vizu(
        stats_data = DAlist_obj(),
        comparison = input$comparison,
        user_gene_list = input$user_gene_list
      )
    })

    # ---------------------- VOLCANO PLOT -----------------------------



    output$input_VP <- plotly::renderPlotly({
      req(stats_plot())
      stats_plot()$volcano_plot
    })

    # ---------------------- PROFILE PLOT -----------------------------

    output$input_PP <- plotly::renderPlotly({
      req(stats_plot(), input$user_gene_list)
      stats_plot()$profile_plot
    })

  })
}

## To be copied in the UI
# mod_glob_stats_explore_ui("glob_stats_explore_1")

## To be copied in the server
# mod_glob_stats_explore_server("glob_stats_explore_1")

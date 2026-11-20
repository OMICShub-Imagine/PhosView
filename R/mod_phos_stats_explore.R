#' phos_stats_explore UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_stats_explore_ui <- function(id) {
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
      column(12, DT::dataTableOutput(ns("PhosTable"))),
      uiOutput(ns("stats_vizu_content"))
    ),

    fluidRow(
      column(4,
             selectizeInput(
               ns("user_phos_list"),
               label = "Select phosphosite to highlight:",
               choices = NULL,
               multiple = FALSE
             ))
    ),

    fluidRow(
      column(6, plotly::plotlyOutput(ns("phos_VP"), height = "400px")),
      column(6, plotly::plotlyOutput(ns("phos_PP"), height = "400px"))
    )

  )
}

#' phos_stats_explore Server Functions
#'
#' @noRd
mod_phos_stats_explore_server <- function(id, DAlist_obj){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    col_info <- reactive({
      req(DAlist_obj())

      if ("Protein" %in% colnames(DAlist_obj()$annotation)) {
        list(
          col_pg   = "Protein",
          col_gene = "Gene.Names",
          col_desc <-  NULL
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
        "user_phos_list",
        choices = sort(DAlist_obj()$annotation$uniprot_id),
        server = TRUE
      )
    })

    # ------------------------- DATATABLE -----------------------------
    output$PhosTable <- DT::renderDataTable({
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

      mod_phos_stats_explore_fct_vizu(
        stats_data = DAlist_obj(),
        comparison = input$comparison,
        user_phos_list = input$user_phos_list
      )
    })

    # ---------------------- VOLCANO PLOT -----------------------------



    output$phos_VP <- plotly::renderPlotly({
      req(stats_plot())
      stats_plot()$volcano_plot
    })

    # ---------------------- PROFILE PLOT -----------------------------

    output$phos_PP <- plotly::renderPlotly({
      req(stats_plot(), input$user_phos_list)
      stats_plot()$profile_plot
    })

  })
}

## To be copied in the UI
# mod_phos_stats_explore_ui("phos_stats_explore_1")

## To be copied in the server
# mod_phos_stats_explore_server("phos_stats_explore_1")

#' merg_ora UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom graphics barplot plot.new text
mod_merg_ora_ui <- function(id) {
  ns <- NS(id)

  ora_tabpanel_factory <- function(ns, title, prefix) {
    tabPanel(

      title = title,

      fluidRow(
        column(12, h3(paste(title, "- UP"))),
        column(12, DT::dataTableOutput(ns(paste0(prefix, "_table_up")))),
        column(12, h3(paste(title, "- DOWN"))),
        column(12, DT::dataTableOutput(ns(paste0(prefix, "_table_down"))))
      ),
      br(),
      fluidRow(
        column(6, plotOutput(ns(paste0(prefix, "_barplot_up")), height = "450px")),
        column(6, plotOutput(ns(paste0(prefix, "_barplot_down")), height = "450px"))
      ),
      br(),
      fluidRow(
        column(6, plotOutput(ns(paste0(prefix, "_cnetplot_up")), height = "450px")),
        column(6, plotOutput(ns(paste0(prefix, "_cnetplot_down")), height = "450px"))
      ),
      br(),
      fluidRow(
        column(6, plotOutput(ns(paste0(prefix, "_emapplot_up")), height = "450px", width = "100%")),
        column(6, plotOutput(ns(paste0(prefix, "_emapplot_down")), height = "450px", width = "100%"))
      )
    )
  }

  tabsetPanel(
    id = ns("ora_tabs"),

    # GLOBAL
    ora_tabpanel_factory(ns, "GO_ORA_Global", "glob_GO"),
    ora_tabpanel_factory(ns, "Reactome_ORA_Global", "glob_Reactome"),
    ora_tabpanel_factory(ns, "KEGG_ORA_Global", "glob_KEGG"),

    # PHOSPHO
    ora_tabpanel_factory(ns, "GO_ORA_Phospho", "phos_GO"),
    ora_tabpanel_factory(ns, "Reactome_ORA_Phospho", "phos_Reactome"),
    ora_tabpanel_factory(ns, "KEGG_ORA_Phospho", "phos_KEGG")
  )



}

#' merg_ora Server Functions
#'
#' @noRd
mod_merg_ora_server <- function(id, ora_results, comparison){
  moduleServer(id, function(input, output, session){

    ###### DATA ACCESS ######

    get_selected <- function(type){
      reactive({
        req(comparison(), ora_results())
        ora_results()[[type]][[comparison()]]
      })
    }

    glob <- get_selected("global")
    phos <- get_selected("phospho")

    ###### GENERIC RENDERERS ######

    render_table <- function(expr) {
      DT::renderDataTable({
        df <- safe_enrich_df(expr())
        DT::datatable(
          df,
          class = "cell-border compact hover order-column",
          rownames = FALSE,
          extensions = c("Scroller", "FixedColumns", "Buttons"),
          selection = "none",
          filter = list(position = 'top'),
          options = list(
            scrollY = 200, scrollX = TRUE, scroller = TRUE,
            fixedColumns = list(leftColumns = 2),
            dom = "frtipB",
            buttons = c('copy','csv','excel','pdf','print','colvis')
          )
        )
      })
    }

    render_barplot <- function(expr) {
      renderPlot({

        req(expr())
        x <- expr()

        if (is.null(x) || is.null(x@result) || nrow(x@result) == 0) {
          plot.new()
          text(0.5, 0.5, "No enrichment available")
          return()
        }

        barplot(x, showCategory = 10)
      })
    }

    render_cnetplot <- function(expr) {
      renderPlot({
        req(expr())

        x <- expr()

        if (is.null(x) ||
            is.null(x@result) ||
            nrow(x@result) == 0) {

          plot.new()
          text(0.5, 0.5, "No enrichment result")

          return()
        }

        # safer threshold: need at least 2 terms
        if (nrow(x@result) < 2) {
          plot.new()
          text(0.5, 0.5, "Not enough enriched terms for network")
          return()
        }

        tryCatch({
          enrichplot::cnetplot(x, showCategory = 10)
        }, error = function(e) {
          plot.new()
          text(0.5, 0.5, "cnetplot failed (no network structure)")
        })
      })
    }

    render_emapplot <- function(expr) {

      renderPlot({

        req(expr())

        x <- expr()

        if (is.null(x) ||
            is.null(x@result) ||
            nrow(x@result) < 2) {

          plot.new()
          text(0.5, 0.5, "Not enough terms for emapplot")
          return()
        }

        sim <- tryCatch(
          enrichplot::pairwise_termsim(x),
          error = function(e) NULL
        )

        if (is.null(sim) ||
            is.null(sim@result) ||
            nrow(sim@result) < 2) {

          plot.new()
          text(0.5, 0.5, "No similarity structure available")
          return()
        }

        tryCatch({

          p <- enrichplot::emapplot(
            sim,
            showCategory = min(10, nrow(sim@result))
          )

          print(p)

        }, error = function(e) {

          message(
            "emapplot failed: ",
            conditionMessage(e)
          )

          plot.new()
          text(
            0.5,
            0.5,
            paste(
              "emapplot failed:",
              conditionMessage(e)
            )
          )
        })

      })
    }

    ###### SAFE BLOCK BUILDER  ######

    build_expr <- function(source, category, direction){

      reactive({
        req(source())

        key <- switch(category,
                      "GO" = paste0("GO_", direction),
                      "Reactome" = paste0("React_", direction),
                      "KEGG" = paste0("KEGG_", direction)
        )

        source()[[key]]
      })
    }

    ###### REGISTER OUTPUTS ######

    register_outputs <- function(source, label){

      get <- function(name){
        reactive({
          req(source())
          source()[[name]]
        })
      }

      ######## TABLES ########

      output[[paste0(label, "_GO_table_up")]]   <- render_table(get("GO_up"))
      output[[paste0(label, "_GO_table_down")]] <- render_table(get("GO_down"))

      output[[paste0(label, "_Reactome_table_up")]]   <- render_table(get("React_up"))
      output[[paste0(label, "_Reactome_table_down")]] <- render_table(get("React_down"))

      output[[paste0(label, "_KEGG_table_up")]]   <- render_table(get("KEGG_up"))
      output[[paste0(label, "_KEGG_table_down")]] <- render_table(get("KEGG_down"))

      ######## BARPLOTS ########

      output[[paste0(label, "_GO_barplot_up")]]   <- render_barplot(get("GO_up"))
      output[[paste0(label, "_GO_barplot_down")]] <- render_barplot(get("GO_down"))

      output[[paste0(label, "_Reactome_barplot_up")]]   <- render_barplot(get("React_up"))
      output[[paste0(label, "_Reactome_barplot_down")]] <- render_barplot(get("React_down"))

      output[[paste0(label, "_KEGG_barplot_up")]]   <- render_barplot(get("KEGG_up"))
      output[[paste0(label, "_KEGG_barplot_down")]] <- render_barplot(get("KEGG_down"))

      ######## CNETPLOT ########

      output[[paste0(label, "_GO_cnetplot_up")]]   <- render_cnetplot(get("GO_up"))
      output[[paste0(label, "_GO_cnetplot_down")]] <- render_cnetplot(get("GO_down"))

      output[[paste0(label, "_Reactome_cnetplot_up")]]   <- render_cnetplot(get("React_up"))
      output[[paste0(label, "_Reactome_cnetplot_down")]] <- render_cnetplot(get("React_down"))

      output[[paste0(label, "_KEGG_cnetplot_up")]]   <- render_cnetplot(get("KEGG_up"))
      output[[paste0(label, "_KEGG_cnetplot_down")]] <- render_cnetplot(get("KEGG_down"))

      ######## EMAPPLOT ########

      output[[paste0(label, "_GO_emapplot_up")]]   <- render_emapplot(get("GO_up"))
      output[[paste0(label, "_GO_emapplot_down")]] <- render_emapplot(get("GO_down"))

      output[[paste0(label, "_Reactome_emapplot_up")]]   <- render_emapplot(get("React_up"))
      output[[paste0(label, "_Reactome_emapplot_down")]] <- render_emapplot(get("React_down"))

      output[[paste0(label, "_KEGG_emapplot_up")]]   <- render_emapplot(get("KEGG_up"))
      output[[paste0(label, "_KEGG_emapplot_down")]] <- render_emapplot(get("KEGG_down"))
    }



    ###### APPLY ######

    register_outputs(glob, "glob")
    register_outputs(phos, "phos")

  })
}


## To be copied in the UI
# mod_merg_ora_ui("merg_ora_1")

## To be copied in the server
# mod_merg_ora_server("merg_ora_1")

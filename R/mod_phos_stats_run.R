#' phos_stats_run UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_stats_run_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # --------------------------------------- #
    # STEP 4 – Action Button to PERFORM STATS #
    # --------------------------------------- #
    fluidRow(
      column(
        12,
        h4("Step 4: Run statistical analysis"),
        actionButton(ns("run_stats"), "PERFORM STATS", class = "btn-lg btn-success"),
        br(), br(),
        verbatimTextOutput(ns("stats_status"))
      )
    )
  )
}

#' phos_stats_run Server Functions
#'
#' @noRd
mod_phos_stats_run_server <- function(id, DAlist_obj, param){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    stats_results <- eventReactive(input$run_stats, {

      validate(
        need(length(param()$contrasts) > 0, "No valid contrasts defined."),
        need(isTruthy(param()$fc_thresh), "FC threshold missing"),
        need(isTruthy(param()$pval_thresh), "p-value threshold missing"),
        need(length(unique(DAlist_obj()$metadata$group)) > 1,
             "Only one group available: cannot build design matrix.")
      )

      tryCatch({

        # ---- Step 1: design ----
        no_intercept <- proteoDA::add_design(
          DAlist_obj(),
          design_formula = ~0 + group
        )

        # ---- Step 2: contrasts ----
        no_intercept <- proteoDA::add_contrasts(
          no_intercept,
          contrasts_vector = param()$contrasts
        )

        # ---- Step 3: fit model ----
        fit <- proteoDA::fit_limma_model(no_intercept)

        # ---- Step 4: extract results ----
        proteoDA::extract_DA_results(
          fit,
          pval_thresh = param()$pval_thresh,
          lfc_thresh = param()$fc_thresh,
          adj_method = "BH",
          extract_intercept = FALSE
        )

      }, error = function(e){
        return(list(error = e$message))
      })

    })

    output$stats_status <- renderPrint({
      out <- stats_results()
      if (is.list(out) && !is.null(out$error)) {
        cat("Error:\n", out$error)
      } else {
        cat("Statistical analysis completed!\n")
      }
    })

    return(stats_results)

  })
}

## To be copied in the UI
# mod_phos_stats_run_ui("phos_stats_run_1")

## To be copied in the server
# mod_phos_stats_run_server("phos_stats_run_1")

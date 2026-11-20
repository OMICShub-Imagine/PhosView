#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  # ---------------------------
  #          APP STATE
  # ---------------------------
  app_state <- reactiveValues(
    mode = "landing",
    session_data = NULL,
    entry_point = NULL
  )

  # ---------------------------
  #        LANDING MODULE
  # ---------------------------
  landing <- mod_landing_server("landing")

  # ---------------------------
  # UI SWITCH
  # ---------------------------
  output$app_ui_main <- renderUI({

    if (app_state$mode == "landing") {

      mod_landing_ui("landing")

    } else {

      # =========================
      #           APP_UI
      # =========================

      navbarPage(
        title = "PhosView",
        id = "phosview",

        tabPanel(
          title = "Global Analysis",
          value = "glob",
          mod_glob_ui("glob")
        ),

        tabPanel(
          title = "Phospho Analysis",
          value = "phos",
          mod_phos_ui("phos")
        ),

        tabPanel(
          title = "Global & Phospho Integration",
          value = "merg",
          mod_merg_ui("merg")
        ),

        tabPanel(
          title = "Export your results",
          value = "export",
          mod_export_ui("export")
        )
      )
    }
  })

  # ---------------------------
  # LANDING ACTIONS
  # ---------------------------

  observeEvent(input$`landing-new_analysis`, {

    app_state$mode <- "analysis"
    app_state$entry_point <- "glob"

  })

  observeEvent(landing$loaded_session(), {

    app_state$session_data <- landing$loaded_session()
    app_state$mode <- "analysis"
    app_state$entry_point <- "merg"

  })

  observeEvent(app_state$mode, {

    req(app_state$mode == "analysis")

    entry <- app_state$entry_point

    session$onFlushed(function() {
      # Disable main tabs
      shinyjs::disable(selector = "a[data-value='glob']")
      shinyjs::disable(selector = "a[data-value='phos']")
      shinyjs::disable(selector = "a[data-value='merg']")
      shinyjs::disable(selector = "a[data-value='export']")

      # Disable all tabs of glob main module except Step 1 at startup
      shinyjs::disable(selector = "a[data-value='glob_filter_sample']")
      shinyjs::disable(selector = "a[data-value='glob_filter_protein']")
      shinyjs::disable(selector = "a[data-value='glob_normalisation']")
      shinyjs::disable(selector = "a[data-value='glob_imputation']")
      shinyjs::disable(selector = "a[data-value='glob_statistics']")

      # Disable tabs of phos main module
      shinyjs::disable(selector = "a[data-value='phos_dalist_creation']")
      shinyjs::disable(selector = "a[data-value='phos_filter_sample']")
      shinyjs::disable(selector = "a[data-value='phos_filter_phos']")
      shinyjs::disable(selector = "a[data-value='phos_normalisation']")
      shinyjs::disable(selector = "a[data-value='phos_imputation']")
      shinyjs::disable(selector = "a[data-value='phos_statistics']")


      updateTabsetPanel(
        session,
        "phosview",
        selected = entry
      )

    })

  })


  glob <- mod_glob_server("glob")

  observeEvent(input[["glob-glob_done"]], {
    shinyjs::enable(selector = "a[data-value='phos']")
    shinyjs::enable(selector = "a[data-value='glob']")
    updateTabsetPanel(session, "phosview", selected = "phos")
  })

  phos <- mod_phos_server("phos", glob = glob$result)

  observeEvent(input[["phos-phos_done"]], {
    shinyjs::enable(selector = "a[data-value='merg']")
    shinyjs::enable(selector = "a[data-value='phos']")
    shinyjs::enable(selector = "a[data-value='glob']")
    updateTabsetPanel(session, "phosview", selected = "merg")
  })

  merg <- mod_merg_server(
    "merg",
    glob = glob,
    phos = phos,
    restored_session = reactive({
      app_state$session_data
    })
  )

  observeEvent(input[["merg-merg_done"]], {
    shinyjs::enable(selector = "a[data-value='export']")
    shinyjs::enable(selector = "a[data-value='merg']")
    shinyjs::enable(selector = "a[data-value='phos']")
    shinyjs::enable(selector = "a[data-value='glob']")
    updateTabsetPanel(session, "phosview", selected = "export")
  })

  mod_export_server(
    "export",
    glob = glob,
    phos = phos,
    ora  = reactive({ merg$ora_results() }),
    kinsub = reactive({ merg$kinsub_results() }),
    pathview = reactive({ merg$pathview_results() }),
    taxo      = merg$taxo
  )

}

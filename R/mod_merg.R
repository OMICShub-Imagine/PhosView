#' merg UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom stats cor na.omit quantile rnorm sd setNames

mod_merg_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        radioButtons(
          inputId = ns("comparison"),
          label = "Select Comparison:",
          choices = c("Waiting for data..." = "placeholder")
        ),

        selectizeInput(
          inputId = ns("user_gene_list"),
          label = "User selected gene(s):",
          choices = NULL,
          selected = NULL,
          multiple = TRUE,
          options = list(maxItems = 1)
        ),

        h4("First perform Over Representation Analysis:"),

        radioButtons(
          inputId = ns("taxo"),
          label = "Choose one of the following taxonomy to perform
          Over Representation Analysis (GO, Reactome and KEGG) on the Global and
          Phospho analysis and click on RUN ORA button:",
          selected = character(0),
          choices = c("human", "mouse"),
          inline = FALSE
        ),

        actionButton(ns("run_ora"), "RUN ORA", class = "btn-lg btn-success"),

        h4("Then perform Kinase-Substrate Analysis:"),

        actionButton(ns("run_signal"), "RUN KinSub", class = "btn-lg btn-success"),

        h4("And finally perform Pathview Analysis:"),

        actionButton(ns("run_pathview"), "RUN PathView", class = "btn-lg btn-success"),

        selectizeInput(
          inputId = ns("kegg_selector"),
          label = "User selected KEGG pathway:",
          choices = NULL,
          selected = NULL,
          multiple = TRUE,
          options = list(maxItems = 1)
        ),

        width = 2
      ),

      mainPanel(
        tabsetPanel(
          id = ns("merg_tabs"),

          tabPanel("DualView",
                   value = "merg_dual_view",
                   mod_merg_dual_view_ui(ns("dual_view"))),

          tabPanel("Over Representation Analysis",
                   value = "merg_ora",
                   mod_merg_ora_ui(ns("ora"))),

          tabPanel("Kinase-Substrate Analysis",
                   value = "merg_signalome",
                   mod_merg_signalome_ui(ns("signalome"))),

          tabPanel("KEGG Pathway Viewer",
                   value = "merg_pathview",
                   mod_merg_pathview_ui(ns("pathview")),
                   fluidRow(
                     column(12, offset = 0, style = "text-align: right;",
                            h4("Move on to the EXPORTS:"),
                            actionButton(ns("merg_done"), "GO to EXPORTS",
                                         class = "btn-lg btn-success")
                     )
                   )
                   )
        )

      )
      )
  )
}

#' merg Server Functions
#'
#' @noRd
mod_merg_server <- function(id, glob, phos, restored_session = reactive(NULL)){
  moduleServer(id, function(input, output, session){

    ns <- session$ns

    #==============================================================
    # SESSION PATHVIEW DIR
    #==============================================================

    pathview_dir <- file.path(
      tempdir(),
      paste0("pathview_", session$token)
    )

    dir.create(
      pathview_dir,
      showWarnings = FALSE,
      recursive = TRUE
    )

    message("PathView session directory: ", pathview_dir)

    #==============================================================
    # DATA ACCESS
    #==============================================================
    glob_data <- reactive({

      sess <- restored_session()

      if (!is.null(sess)) {
        return(sess$glob$stats_data)
      }

      req(glob$result())
      glob$result()$stats_data
    })

    phos_data <- reactive({

      sess <- restored_session()

      if (!is.null(sess)) {
        return(sess$phos$stats_data)
      }

      req(phos$result())
      phos$result()$stats_data
    })

    #==============================================================
    # UI INIT
    #==============================================================
    observe({

      req(glob_data(), phos_data())

      results <- glob_data()$results

      if (length(results) > 0) {
        updateRadioButtons(
          session,
          inputId = "comparison",
          choices = names(results),
          selected = names(results)[1]
        )
      }

      # Use annotation$gene (already harmonized upstream)
      updateSelectizeInput(
        session,
        "user_gene_list",
        choices = sort(unique(glob_data()$annotation$gene)),
        server = TRUE
      )
    })

    #==============================================================
    # CALL MODULES ONCE (NOT INSIDE OBSERVERS)
    #==============================================================
    mod_merg_dual_view_server(
      "dual_view",
      glob = glob,
      phos = phos,
      comparison = reactive(input$comparison),
      user_gene_list = reactive(input$user_gene_list),
      restored_session = restored_session
    )

    #==============================================================
    # ORA
    #==============================================================
    ora_results <- reactiveVal(NULL)

    observe({

      sess <- restored_session()

      req(sess)

      if (!is.null(sess$ora)) {
        ora_results(sess$ora)
      }

    })

    observeEvent(input$run_ora, {

      req(input$taxo, glob_data(), phos_data())

      if (input$taxo == "human") {
        AnnotDb <- org.Hs.eg.db::org.Hs.eg.db
        taxo_reactome <- "human"
        taxo_kegg <- "hsa"
      } else {
        AnnotDb <- org.Mm.eg.db::org.Mm.eg.db
        taxo_reactome <- "mouse"
        taxo_kegg <- "mmu"
      }

      stats_global  <- glob_data()
      stats_phospho <- phos_data()

      withProgress(message = "Performing ORA analysis...", value = 0, {

        incProgress(0.1, "Running global ORA...")

        ora_global <- run_all_ORA(
          stats_data = stats_global,
          AnnotDb = AnnotDb,
          taxo_reactome = taxo_reactome,
          taxo_kegg = taxo_kegg,
          mode = "global"
        )

        incProgress(0.5, "Running phospho ORA...")

        ora_phospho <- run_all_ORA(
          stats_data = stats_phospho,
          AnnotDb = AnnotDb,
          taxo_reactome = taxo_reactome,
          taxo_kegg = taxo_kegg,
          mode = "phospho"
        )

        incProgress(1, "Done")
      })

      ora_results(list(
        global = ora_global,
        phospho = ora_phospho
      ))
    })

    # Call mod_merg_ora_serv
    mod_merg_ora_server(
      "ora",
      ora_results = ora_results,
      comparison = reactive(input$comparison)
    )

    #==============================================================
    # KINASE-SUBSTRATE
    #==============================================================
    kinsub_results <- reactiveVal(NULL)

    observe({

      sess <- restored_session()

      req(sess)

      if (!is.null(sess$kinsub)) {
        kinsub_results(sess$kinsub)
      }

    })

    observeEvent(input$run_signal, {

      req(input$taxo, phos_data())

      if (input$taxo == "human") {
        species <- "human"
        PhosphoSite.Db <- PhosphoSite.human
      } else {
        species <- "mouse"
        PhosphoSite.Db <- PhosphoSite.mouse
      }

      stats_phospho <- phos_data()

      withProgress(message = "Kinase-Substrate analysis...", value = 0, {

        incProgress(0.5)

        kinsub_res <- compute_signalome(
          phospho_stats = stats_phospho,
          species = species,
          PhosphoSite.Db = PhosphoSite.Db
        )

        incProgress(1)
      })

      kinsub_results(kinsub_res)

      message("=== KINSUB STORED ===")
      print(str(kinsub_results()))
    })

    mod_merg_signalome_server(
      "signalome",
      kinsub_results = kinsub_results,
      comparison = reactive(input$comparison)
    )

    #####################################
    ######## PATHVIEW ANALYSIS ##########
    #####################################

    pathview_results <- reactiveVal(NULL)
    saved_taxo <- reactiveVal(NULL)
    taxo_restore_pending <- reactiveVal(FALSE)

    # ---- RESTORE SAVED PATHVIEW RESULTS ----

    # ---- RESTORE SAVED PATHVIEW RESULTS ----

    observeEvent(restored_session(), {

      message("\n=== RESTORATION OBSERVER EXECUTED ===")

      sess <- restored_session()
      req(sess)

      message("Session fields: ", paste(names(sess), collapse = ", "))
      message("Saved taxo: '", sess$taxo, "'")

      message("=== RESTORING PATHVIEW SESSION ===")

      # 1. Restore PathView results first
      if (!is.null(sess$pathview)) {
        pathview_results(sess$pathview)
        message("Saved PathView results restored")
      }

      # 2. Recover taxonomy from saved session
      restored_taxo <- sess$taxo

      # The saved taxonomy may be NULL or an empty string
      if (is.null(restored_taxo) ||
          length(restored_taxo) == 0 ||
          is.na(restored_taxo[1]) ||
          trimws(restored_taxo[1]) == "") {

        message("Saved taxonomy is empty. Inferring from KEGG results...")

        restored_taxo <- NULL

        # Inspect saved KEGG enrichment results
        if (!is.null(sess$pathview$enrichment)) {

          kegg_ids <- unlist(
            lapply(sess$pathview$enrichment, function(x) {

              if (!is.null(x$KEGG) &&
                  !is.null(x$KEGG@result) &&
                  nrow(x$KEGG@result) > 0) {

                return(x$KEGG@result$ID)
              }

              return(character(0))
            }),
            use.names = FALSE
          )

          message("KEGG pathway ID examples: ",
                  paste(head(kegg_ids, 10), collapse = ", "))

          if (any(grepl("^hsa", kegg_ids))) {

            restored_taxo <- "human"

          } else if (any(grepl("^mmu", kegg_ids))) {

            restored_taxo <- "mouse"
          }
        }
      }

      message("Recovered taxonomy: ", restored_taxo)

      # 3. Restore taxonomy radio button
      if (!is.null(restored_taxo) &&
          length(restored_taxo) == 1 &&
          !is.na(restored_taxo) &&
          restored_taxo %in% c("human", "mouse")) {

        saved_taxo(restored_taxo)

        taxo_restore_pending(
          !identical(input$taxo, restored_taxo)
        )

        updateRadioButtons(
          session,
          "taxo",
          selected = restored_taxo
        )

        message("Taxonomy restoration requested: ", restored_taxo)

      } else {

        taxo_restore_pending(FALSE)

        warning(
          "Could not recover taxonomy from the saved session. ",
          "Please check the saved KEGG enrichment results."
        )
      }

      # 4. Restore comparison
      if (!is.null(sess$comparison)) {

        updateRadioButtons(
          session,
          "comparison",
          selected = sess$comparison
        )
      }

      message("=== PATHVIEW SESSION RESTORATION COMPLETED ===")

    }, ignoreInit = FALSE, once = TRUE)


    observeEvent(input$taxo, {
      message("=== TAXONOMY INPUT CHANGED ===")
      message("New input$taxo: '", input$taxo, "'")
    }, ignoreInit = TRUE)


    # ---- WAIT FOR TAXONOMY RESTORATION ----

    observeEvent(input$taxo, {

      if (isTRUE(taxo_restore_pending()) &&
          identical(input$taxo, saved_taxo())) {

        taxo_restore_pending(FALSE)

        message("Taxonomy restoration completed")
      }

    }, ignoreInit = TRUE)


    observeEvent(input$run_pathview, {

      #####################################################
      ###### FIRST STEP: RUN THE PATHVIEW ANALYSIS ########
      #####################################################

      message("=== PATHVIEW BUTTON CLICKED ===")

      req(input$taxo, phos_data(), kinsub_results())

      message("Inputs validated")

      if (input$taxo == "human") {
        AnnotDb <- org.Hs.eg.db::org.Hs.eg.db
        taxo_kegg <- "hsa"
      } else {
        AnnotDb <- org.Mm.eg.db::org.Mm.eg.db
        taxo_kegg <- "mmu"
      }

      message("Species selected: ", taxo_kegg)

      stats_phospho   <- phos_data()
      kinsub_res_data <- kinsub_results()

      message("Data extracted")

      withProgress(message = "Pathview analysis...", value = 0, {

        message("Before compute_pathview()")

        incProgress(0.5)

        pathview_res <- compute_pathview(
          stat_data   = stats_phospho,
          kinsub_data = kinsub_res_data,
          AnnotDb     = AnnotDb,
          taxo_kegg   = taxo_kegg
        )

        # NEW DEBUG LAYER
        message("=== SHINY OUTPUT CHECK ===")
        str(pathview_res)
        stopifnot(!is.null(pathview_res$final_input))
        stopifnot(nrow(pathview_res$final_input) > 0)

        message("After compute_pathview()")

        incProgress(1)
      })

      message("Before pathview_results() assignment")

      pathview_results(pathview_res)

      message("After pathview_results() assignment")
    })

    ###############################################################
    ########## SECOND STEP: DISPLAY KEGG PATHWAY ##################
    ###############################################################

    # Selected comparison ORA results
    kinsub_ora_selected_results <- reactive({
      req(pathview_results(), input$comparison)
      pathview_results()$enrichment[[input$comparison]]
    })

    kinsub_selected_results <- reactive({
      req(pathview_results(), input$comparison)
      pathview_results()$final_input %>%
        dplyr::filter(Comparison %in% input$comparison)
    })

    mod_merg_pathview_server(
      "pathview",

      kinsub_ora_selected_results = kinsub_ora_selected_results,
      kinsub_selected_results     = kinsub_selected_results,
      comparison                  = reactive(input$comparison),
      kegg_selector               = reactive(input$kegg_selector),


      taxo_kegg = reactive({

        req(!taxo_restore_pending())
        req(input$taxo)

        if (input$taxo == "human") {
          "hsa"
        } else {
          "mmu"
        }

      }),

      AnnotDb = reactive({

        message("=== DEBUG AnnotDb reactive ===")
        message("input$taxo: ", input$taxo)
        message("taxo_restore_pending: ", taxo_restore_pending())

        req(!taxo_restore_pending())
        req(input$taxo)

        if (input$taxo == "human") {
          message("Returning human annotation database")
          org.Hs.eg.db::org.Hs.eg.db
        } else {
          message("Returning mouse annotation database")
          org.Mm.eg.db::org.Mm.eg.db
        }
      }),


      pathview_dir = pathview_dir
    )

    ###############################################################
    ################# UPDATE KEGG SELECTOR  #######################
    ###############################################################

    #==============================================================
    # UPDATE KEGG SELECTOR AFTER SESSION RESTORATION
    #==============================================================

    observeEvent(
      list(pathview_results(), input$comparison),
      {

        req(pathview_results(), input$comparison)

        message("=== UPDATING KEGG SELECTOR ===")
        message("Comparison: ", input$comparison)

        comp_res <- pathview_results()$enrichment[[input$comparison]]

        if (is.null(comp_res) ||
            is.null(comp_res$KEGG) ||
            is.null(comp_res$KEGG@result) ||
            nrow(comp_res$KEGG@result) == 0) {

          updateSelectizeInput(
            session,
            "kegg_selector",
            choices = character(0),
            selected = character(0)
          )

          message("No KEGG enrichment available")

          return()
        }

        kegg_res <- comp_res$KEGG@result

        kegg_choices <- setNames(
          kegg_res$ID,
          kegg_res$Description
        )

        message("Available KEGG pathways: ", length(kegg_choices))

        # Select the first available pathway
        selected_pathway <- kegg_res$ID[1]

        updateSelectizeInput(
          session,
          "kegg_selector",
          choices = kegg_choices,
          selected = selected_pathway,
          server = TRUE
        )

        message("Selected KEGG pathway: ", selected_pathway)

      },

      ignoreInit = FALSE
    )


    #==============================================================
    # EXPORT
    #==============================================================
    observeEvent(input$merg_done, {
      session$sendCustomMessage("merg_finished", TRUE)
    })

    return(list(
      ora_results      = ora_results,
      kinsub_results   = kinsub_results,
      pathview_results = pathview_results
    ))

  })
}

## To be copied in the UI
# mod_merg_ui("merg_1")

## To be copied in the server
# mod_merg_server("merg_1")

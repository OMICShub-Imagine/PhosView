#' phos_stats_param UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_stats_param_ui <- function(id) {
  ns <- NS(id)
  tagList(

    # ------------------------ #
    # STEP 1 & 2 – Comparisons #
    # ------------------------ #
    fluidRow(
      column(
        12,
        h4("Step 1: Set the number of statistical comparisons:"),
        numericInput(ns("n"), "Number of comparisons", value = 1, min = 1),
        h4("Step 2: Fill the contrast information:"),
        uiOutput(ns("comparisons_ui")),
        checkboxInput(ns("auto_name"), "Auto-generate contrast names", TRUE)
      )
    ),

    # ----------------------------- #
    # STEP 3 – Set Stats Thresholds #
    # ----------------------------- #
    fluidRow(
      column(
        12,
        h4("Step 3: Set statistical thresholds")
      ),
      column(
        6,
        numericInput(ns("fc_thresh"), "Fold-change threshold (log2FC)", value = 1, step = 0.1)
      ),
      column(
        6,
        numericInput(ns("pval_thresh"), "Adjusted p-value threshold", value = 0.05, step = 0.01)
      )
    )

  )
}

#' phos_stats_param Server Functions
#'
#' @noRd
mod_phos_stats_param_server <- function(id, DAlist_obj){

  moduleServer(id, function(input, output, session){
    ns <- session$ns

    # ------ #
    # GROUPS #
    # ------ #
    groups <- reactive({
      req(DAlist_obj())
      unique(DAlist_obj()$metadata$group)
    })

    # ---------------------------- #
    # STEP 1 – Comparison/Contrast #
    # ---------------------------- #
    output$comparisons_ui <- renderUI({
      req(input$n, groups())

      tagList(
        lapply(seq_len(input$n), function(i){
          fluidRow(
            column(
              4,
              textInput(
                ns(paste0("comp_name_", i)),
                paste("Comparison name", i),
                placeholder = "e.g. KO_vs_WT"
              )
            ),
            column(
              4,
              selectInput(
                ns(paste0("left_", i)),
                paste("Left group", i),
                choices = groups()
              )
            ),
            column(
              4,
              selectInput(
                ns(paste0("right_", i)),
                paste("Right group", i),
                choices = groups()
              )
            )
          )
        })
      )
    })

    # -------------------------------------------- #
    # AUTO-GENERATE CONTRAST NAMES (SAFE VERSION) #
    # -------------------------------------------- #
    observe({
      req(input$n)

      lapply(seq_len(input$n), function(i){

        left  <- input[[paste0("left_", i)]]
        right <- input[[paste0("right_", i)]]

        # must have both selected
        if (!isTruthy(left) || !isTruthy(right)) return()
        if (left == right) return()

        new_name <- paste0(left, "_vs_", right)

        current <- input[[paste0("comp_name_", i)]]

        # only update if:
        # - empty OR
        # - already auto-generated OR
        # - matches reversed version
        if (is.null(current) ||
            current == "" ||
            grepl("_vs_", current)) {

          updateTextInput(
            session,
            paste0("comp_name_", i),
            value = new_name
          )
        }

      })
    })

    # ---------------------------- #
    # CONTRASTS VALIDATION + BUILD #
    # ---------------------------- #
    contrasts_reactive <- reactive({
      req(input$n)

      contrasts <- character(input$n)
      warnings  <- character(0)

      for(i in seq_len(input$n)){
        comp  <- input[[paste0("comp_name_", i)]]
        left  <- input[[paste0("left_", i)]]
        right <- input[[paste0("right_", i)]]

        if (isTruthy(comp) && isTruthy(left) && isTruthy(right)) {

          if (left == right) {
            warnings <- c(warnings, paste0("Contrast ", comp, " ignored (left == right)"))
            next
          }

          contrasts[i] <- paste0(comp, " = ", left, " - ", right)
        }
      }

      if (length(warnings) > 0)
        message(paste(warnings, collapse = "\n"))

      contrasts[contrasts != ""]
    })

    # ---------------------------- #
    # RETURN ALL PARAMETERS         #
    # ---------------------------- #
    return(
      reactive(list(
        contrasts   = contrasts_reactive(),
        fc_thresh   = input$fc_thresh,
        pval_thresh = input$pval_thresh
      ))
    )

  })
}

## To be copied in the UI
# mod_phos_stats_param_ui("phos_stats_param_1")

## To be copied in the server
# mod_phos_stats_param_server("phos_stats_param_1")

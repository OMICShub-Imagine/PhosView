#' phos_first_class_filt UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_first_class_filt_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             h4("Step 4: According to the QC report generated just above,
                 how would you like to filter your dataset (row-wise) based on the presence
                 of phosphosite having localisation probaiblity >75% in at least 50% or 70% or 100% of replicates
                 in at least one of your experimental group"),

             radioButtons(
               inputId = ns("remove_phos"),
               label = "Do you want to keep phosphosites that have localisation probability >75% in at least:",
               selected = character(0),
               choices = c("50% of replicates in at least one of your experimental group",
                           "70% of replicates in at least one of your experimental group",
                           "100% of replicates in at least one of your experimental group"),
               inline = TRUE
             ),

             uiOutput(ns("keep_button")),

             verbatimTextOutput(ns("upload_status"))
      )
    )
  )
}

#' phos_first_class_filt Server Functions
#'
#' @noRd
mod_phos_first_class_filt_server <- function(id, quant_data, meta_data){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    # Reactive value to hold filtered object
    first_class_phos <- reactiveVal(NULL)


    # Dynamically render REMOVE button
    output$keep_button <- renderUI({
      req(input$remove_phos)
      actionButton(ns("keep"), "KEEP", class = "btn-lg btn-success")
    })

    # Observe remove button
    observeEvent(input$keep, {
      threshold_id <- input$remove_phos
      if (is.null(threshold_id) || length(threshold_id) == 0) {
        showNotification("No threshold selected!", type = "warning")
        return()
      }

      tryCatch({

        # Filter using your helper function
        filt_phos <- mod_phos_first_class_filt_fct_transform(
          quant_data = quant_data(),
          meta_data = meta_data(),
          threshold = threshold_id
        )

        first_class_phos(filt_phos)
        output$upload_status <- renderText("PTMsite matrix successfully filtered row-wise.")

        # Informative message
        removed_n <- nrow(quant_data()) - nrow(filt_phos)
        kept_n <- nrow(filt_phos)
        showNotification(
          paste0("Removed ", removed_n, " PTMsites; ",
                 kept_n, " first class phosphosite kept (", round(kept_n / (kept_n + removed_n) * 100, 1), "% retained)."),
          type = "message"
        )

      },
      error = function(e) {
        msg <- paste("Error during filtering:", e$message)
        output$upload_status <- renderText(msg)
        showNotification(msg, type = "error")
      })
    })


    # Return the filtered object
    return(first_class_phos)

  })
}

## To be copied in the UI
# mod_phos_first_class_filt_ui("phos_first_class_filt_1")

## To be copied in the server
# mod_phos_first_class_filt_server("phos_first_class_filt_1")

#' phos_filt_phos_transform UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_phos_filt_phos_transform_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             h4("Step 1: According to the QC report generated in the Filter sample(s) step,
               how would you like to filter your dataset (row-wise) based on the presence
               of 50% or 70% or 100% valid values (VV) in at least one of your experimental group"),

             radioButtons(
               inputId = ns("remove_phos"),
               label = "Do you want to keep phosphosites that have at least:",
               selected = character(0),
               choices = c("50% of VV in at least one of your experimental group",
                           "70% of VV in at least one of your experimental group",
                           "100% of VV in at least one of your experimental group"),
               inline = TRUE
             ),

             uiOutput(ns("keep_button")),

             verbatimTextOutput(ns("upload_status"))
      )
    )

  )
}

#' phos_filt_phos_transform Server Functions
#'
#' @noRd
mod_phos_filt_phos_transform_server <- function(id, DAlist_obj){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    # Reactive value to hold filtered object
    filter_phos_obj <- reactiveVal(NULL)


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
        req(DAlist_obj())
        # Filter using your helper function
        filt_phos <- mod_phos_filt_phos_transform_fct_filt_phos(
          DAList_obj = DAlist_obj(),
          threshold = threshold_id
        )

        filter_phos_obj(filt_phos)
        output$upload_status <- renderText("DAlist successfully filtered row-wise.")

        # Informative message
        removed_n <- nrow(DAlist_obj()$data) - nrow(filt_phos$data)
        kept_n <- nrow(filt_phos$data)
        showNotification(
          paste0("Removed ", removed_n, " phosphosites; ",
                 kept_n, " phosphosites kept (", round(kept_n / (kept_n + removed_n) * 100, 1), "% retained)."),
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
    return(filter_phos_obj)

  })
}

## To be copied in the UI
# mod_phos_filt_phos_transform_ui("phos_filt_phos_transform_1")

## To be copied in the server
# mod_phos_filt_phos_transform_server("phos_filt_phos_transform_1")

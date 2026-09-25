#' glob_filt_prot_transform UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_glob_filt_prot_transform_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
    column(12,
           h4("Step 1: According to the QC report generated in the Filter sample(s) step,
               how would you like to filter your dataset (row-wise) based on the presence
               of 50% or 70% or 100% valid values (VV) in at least one of your experimental group"),

           radioButtons(
             inputId = ns("remove_prot"),
             label = "Do you want to keep proteins that have at least:",
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

#' glob_filt_prot_transform Server Functions
#'
#' @noRd
mod_glob_filt_prot_transform_server <- function(id, DAlist_obj){
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive value to hold filtered object
    filter_prot_obj <- reactiveVal(NULL)


    # Dynamically render REMOVE button
    output$keep_button <- renderUI({
      req(input$remove_prot)
      actionButton(ns("keep"), "KEEP", class = "btn-lg btn-success")
    })

    # Observe remove button
    observeEvent(input$keep, {
      threshold_id <- input$remove_prot
      if (is.null(threshold_id) || length(threshold_id) == 0) {
        showNotification("No threshold selected!", type = "warning")
        return()
      }

      tryCatch({
        req(DAlist_obj())
        # Filter using your helper function
        filt_prot <- mod_glob_filt_prot_transform_fct_filt_prot(
          DAList_obj = DAlist_obj(),
          threshold = threshold_id
        )

        filter_prot_obj(filt_prot)
        output$upload_status <- renderText("DAlist successfully filtered row-wise.")

        # Informative message
        removed_n <- nrow(DAlist_obj()$data) - nrow(filt_prot$data)
        kept_n <- nrow(filt_prot$data)
        showNotification(
          paste0("Removed ", removed_n, " proteins; ",
                 kept_n, " proteins kept (", round(kept_n / (kept_n + removed_n) * 100, 1), "% retained)."),
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
    return(filter_prot_obj)
  })
}

## To be copied in the UI
# mod_glob_filt_prot_transform_ui("glob_filt_prot_transform_1")

## To be copied in the server
# mod_glob_filt_prot_transform_server("glob_filt_prot_transform_1")

#' glob_filt_sample_transform UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList

mod_glob_filt_sample_transform_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             h4("Step 1: According to the QC report generated in the DAlist creation step,
                you have the opportunity to remove outlier(s) and therefore filter your dataset column-wise"),

             radioButtons(
               inputId = ns("remove_outlier"),
               label = "Do you want to remove any outlier?",
               selected = character(0),
               choices = c("No", "Yes"),
               inline = TRUE
             ),

             h4("Step 2: If you click on Yes, choose any outlier(s) in the following list. If you click on No, go directly to Step3"),
             uiOutput(ns("outlier_select_ui")),

             uiOutput(ns("remove_button")),

             verbatimTextOutput(ns("upload_status"))
      )
    )
  )
}


#' glob_filt_sample_transform Server Function
#'
#' @noRd
mod_glob_filt_sample_transform_server <- function(id, DAlist_obj) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive value to hold filtered object
    filt_sample_obj <- reactiveVal(NULL)

    # Extract sample names
    sample_names <- reactive({
      req(DAlist_obj())
      DAlist_obj()$metadata$Simplify_ID
    })

    # Dynamically render selectInput
    output$outlier_select_ui <- renderUI({
      req(input$remove_outlier)
      if (input$remove_outlier == "Yes") {
        selectInput(
          inputId = ns("outlier_samples"),
          label = "Select sample(s) to remove:",
          choices = sample_names(),
          multiple = TRUE,
          selectize = TRUE
        )
      }
    })

    # Dynamically render REMOVE button
    output$remove_button <- renderUI({
      req(input$remove_outlier)
      if (input$remove_outlier == "Yes") {
        actionButton(ns("remove"), "REMOVE", class = "btn-lg btn-success")
      }
    })

    # Observe remove button
    observeEvent(input$remove, {
      outliers_id <- input$outlier_samples
      if (is.null(outliers_id) || length(outliers_id) == 0) {
        showNotification("No sample selected!", type = "warning")
        return()
      }

      tryCatch({
        req(DAlist_obj())
        # Filter using your dedicated helper function
        filt_samples <- mod_glob_filt_sample_transform_fct_filt_samples(
          DAlist_obj(),
          outliers_id
        )

        filt_sample_obj(filt_samples)
        output$upload_status <- renderText("DAlist successfully filtered column-wise.")
        showNotification(paste("Removed:", paste(outliers_id, collapse = ", ")), type = "message")
      },
      error = function(e) {
        msg <- paste("Error during filtering:", e$message)
        output$upload_status <- renderText(msg)
        showNotification(msg, type = "error")
      })
    })

    # If user selects "No", keep original object
    observe({
      req(input$remove_outlier)
      if (input$remove_outlier == "No") {
        filt_sample_obj(DAlist_obj())
      }
    })

    # Return the filtered object
    return(filt_sample_obj)
  })
}



## To be copied in the UI
# mod_glob_filt_sample_transform_ui("glob_filt_sample_transform_1")

## To be copied in the server
# mod_glob_filt_sample_transform_server("glob_filt_sample_transform_1")

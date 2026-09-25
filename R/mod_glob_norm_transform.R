#' glob_norm_transform UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_glob_norm_transform_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             h4("Step 2: According to the normalization report display above, would you like to apply
                any normalization (followed by log2 transformation of your dataset) or just log2 transform your dataset"),

             radioButtons(
               inputId = ns("apply_norm"),
               label = "Do you want to normalize your data?",
               selected = character(0),
               choices = c("No", "Yes"),
               inline = TRUE
             ),

             h4("Step 3: If you click on Yes, choose one of the following normalization methods."),
             uiOutput(ns("norm_select_ui")),

             uiOutput(ns("norm_button")),

             verbatimTextOutput(ns("upload_status"))
      )
    )

  )
}

#' glob_norm_transform Server Functions
#'
#' @noRd
mod_glob_norm_transform_server <- function(id, DAlist_obj){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    # Reactive value to hold filtered object
    norm_data_obj <- reactiveVal(NULL)

    # Dynamically render selectInput
    output$norm_select_ui <- renderUI({
      req(input$apply_norm)
      if (input$apply_norm == "Yes") {
        radioButtons(
          inputId = ns("norm_method"),
          label = "Select one of the following normalization method:",
          choices = c("median normalization followed by log2 transformation",
                      "mean normalization followed by log2 transformation",
                      "vsn (variance stabilizing normalization) normalization followed by log2 transformation",
                      "quantile normalization followed by log2 transformation",
                      "cycloess (cyclic loess) normalization followed by log2 transformation",
                      "rlr (robust linear regression) normalization followed by log2 transformation",
                      "gi (global intensity) normalization followed by log2 transformation"),
          selected = character(0),
          inline = TRUE
        )
      }
      else {
        radioButtons(
          inputId = ns("norm_method"),
          label = "Apply only log2 transformation",
          choices = c("no normalization but log2 transformation"),
          selected = character(0),
          inline = TRUE
        )
          }
    })

    # Dynamically render REMOVE button
    output$norm_button <- renderUI({
      req(input$apply_norm)
      actionButton(ns("norm"), "NORMALIZATION", class = "btn-lg btn-success")
    })

    # Observe remove button
    observeEvent(input$norm, {
      norm_id <- input$norm_method
      if (is.null(norm_id) || length(norm_id) == 0) {
        showNotification("No normalization method selected!", type = "warning")
        return()
      }

      tryCatch({
        req(DAlist_obj())
        # Filter using your dedicated helper function
        norm_data <- mod_glob_norm_transform_fct_norm(
          DAList_obj = DAlist_obj(),
          norm = norm_id
        )

        norm_data_obj(norm_data)
        output$upload_status <- renderText("DAlist successfully normalized!")
        showNotification(paste("Apply:", paste(norm_id)), type = "message")
      },
      error = function(e) {
        msg <- paste("Error during normalization:", e$message)
        output$upload_status <- renderText(msg)
        showNotification(msg, type = "error")
      })
    })


    # Return the filtered object
    return(norm_data_obj)

  })
}

## To be copied in the UI
# mod_glob_norm_transform_ui("glob_norm_transform_1")

## To be copied in the server
# mod_glob_norm_transform_server("glob_norm_transform_1")

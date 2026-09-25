#' landing UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_landing_ui <- function(id){

  ns <- NS(id)

  tagList(

    # ---- HERO SECTION ----
    div(
      style = "
        min-height: 100vh;
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
        background: linear-gradient(135deg, #f8fafc, #eef2ff);
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Arial;
      ",

      div(
        style = "
          width: 1000x;
          padding: 32px;
          border-radius: 16px;
          background: white;
          box-shadow: 0 10px 30px rgba(0,0,0,0.12);
          text-align: center;
        ",

        # ---- LOGO ----
        tags$img(
          src = "www/PhosView_landing_page.png",
          style = "width: 60%;
                   max-width: 300px;
                   height: auto;
                   margin-bottom: 20px;"
        ),

        # ---- DESCRIPTION ----
        div(
          "PhosView is a Shiny app that allows you to:",
          tags$ul(
            tags$li("process your global proteomics data (i.e. Global Analysis tab)"),
            tags$li("process your phosphoproteomics data (i.e. Phospho Analysis tab)"),
            tags$li("integrate your global- and phospho- proteomics data together (i.e. Global & Phospho Analysis tab)"),
            tags$li("export your results (i.e. Export your results tab)")
          ),
          style = "font-size: 15px; color: #666; margin-bottom: 25px;text-align: left;"
        ),

        # ---- CTA BUTTON ----
        actionButton(
          ns("new_analysis"),
          "New Analysis",
          style = "
            width: 100%;
            background: #198754;
            color: white;
            border: none;
            padding: 12px;
            font-size: 16px;
            border-radius: 8px;
            margin-bottom: 15px;
          "
        ),

        # ---- FILE INPUT ----
        fileInput(
          ns("session_file"),
          "Reload Session",
          width = "100%"
        ),

        br(),

        actionButton(
          ns("load_session"),
          "Load Session",
          style = "
            width: 100%;
            background: #198754;
            color: white;
            border: none;
            padding: 12px;
            font-size: 16px;
            border-radius: 8px;
            margin-bottom: 15px;
          "
        )
      )
    )
  )
}

#' landing Server Functions
#'
#' @noRd
mod_landing_server <- function(id){

  moduleServer(id, function(input, output, session){

    loaded_session <- reactiveVal(NULL)

    observeEvent(input$load_session, {

      req(input$session_file)

      obj <- readRDS(input$session_file$datapath)

      loaded_session(obj)

    })

    return(list(
      loaded_session = loaded_session
    ))

  })

}

## To be copied in the UI
# mod_landing_ui("landing_1")

## To be copied in the server
# mod_landing_server("landing_1")

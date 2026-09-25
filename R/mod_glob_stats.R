#' glob_stats UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_glob_stats_ui <- function(id) {
  ns <- NS(id)
  tagList(
    mod_glob_stats_param_ui(ns("param")),
    mod_glob_stats_run_ui(ns("run")),
    mod_glob_stats_explore_ui(ns("explore"))

  )
}

#' glob_stats Server Functions
#'
#' @noRd
mod_glob_stats_server <- function(id, imput_data){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    # Step 1: User-defined contrasts + thresholds
    stats_param <- mod_glob_stats_param_server(
      id = "param",
      DAlist_obj = imput_data
    )

    # Step 2: Run statistical analysis
    stats_data <- mod_glob_stats_run_server(
      id = "run",
      DAlist_obj = imput_data,
      param = stats_param
    )

    # Step 3: Explore results (volcano, profile, table)
    stats_explore <- mod_glob_stats_explore_server(
      id = "explore",
      DAlist_obj = stats_data
    )

    return(stats_data)

  })
}

## To be copied in the UI
# mod_glob_stats_ui("glob_stats_1")

## To be copied in the server
# mod_glob_stats_server("glob_stats_1")

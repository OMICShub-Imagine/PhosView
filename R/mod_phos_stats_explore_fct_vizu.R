#' Generate volcano and profile plots for phosphosite statistics
#'
#' This function generates interactive visualizations of phosphosite statistics for a given comparison:
#' a volcano plot highlighting up-/down-regulated sites, and an optional profile plot for selected phosphosites.
#'
#' @param stats_data A list returned by your statistical analysis module. It should contain:
#'   \itemize{
#'     \item \code{results}: named list of data frames with differential analysis results per comparison
#'     \item \code{data}: matrix or data frame of LFQ intensity values
#'     \item \code{annotation}: data frame mapping phosphosites to genes and protein descriptions
#'     \item \code{metadata}: data frame with sample information (e.g., group, Simplify_ID)
#'   }
#' @param comparison Character string specifying which comparison to visualize (must match a name in \code{stats_data$results})
#' @param user_phos_list Optional character vector of phosphosites (uniprot_id) to highlight in the volcano plot and to generate profile plots.
#'                        If NULL or empty, no profile plot is generated.
#'
#' @return A named list containing:
#' \describe{
#'   \item{volcano_plot}{Interactive \pkg{plotly} volcano plot highlighting up-/down-regulated and optionally selected phosphosites}
#'   \item{profile_plot}{Interactive \pkg{plotly} intensity profile plot for the selected phosphosite(s). NULL if no phosphosite is selected}
#' }
#'
#' @import ggplot2
#' @import plotly
#' @importFrom dplyr left_join case_when
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @examples
#' \dontrun{
#' plots <- mod_phos_stats_explore_fct_vizu(stats_data = my_stats, comparison = "Treatment_vs_Control", user_phos_list = "P12345")
#' plotly::ggplotly(plots$volcano_plot)
#' plotly::ggplotly(plots$profile_plot)
#' }

mod_phos_stats_explore_fct_vizu <- function(stats_data, comparison, user_phos_list = NULL) {

  if ("Protein" %in% colnames(stats_data$annotation)) {
    col_pg <- "Protein"
    col_gene <- "Gene.Names"
    col_desc <-  NULL
  } else if ("PG.ProteinGroups" %in% colnames(stats_data$annotation)) {
    col_pg <- "PG.ProteinGroups"
    col_gene <- "PG.Genes"
    col_desc <- "PG.ProteinDescriptions"
  }

  ############################ VOLCANO PLOT #########################################

  df_volcano <- stats_data$results[[comparison]] %>%
    tibble::rownames_to_column("uniprot_id") %>%
    dplyr::left_join(., stats_data$annotation, by = "uniprot_id")

  req(df_volcano)

  # Assign significance categories
  df_volcano$point_color <- dplyr::case_when(
    df_volcano$sig.FDR == -1 ~ "DownRegulated",
    df_volcano$sig.FDR == 1  ~ "UpRegulated",
    TRUE                     ~ "nonSignificant"
  )

  # Highlight selected gene if provided
  if (!is.null(user_phos_list) & user_phos_list != "") {
    df_volcano$point_color[df_volcano$uniprot_id == user_phos_list] <- "selected"
  }

  colors_volcano <- c(
    "DownRegulated" = "#5555FF",
    "nonSignificant" = "grey",
    "UpRegulated" = "#FF5555",
    "selected" = "black"
  )

  # Build tooltip dynamically
  tooltip_text <- paste0(
    "Phosphosite: ", df_volcano$uniprot_id,
    "<br>Gene: ", df_volcano[[col_gene]]
  )

  if (!is.null(col_desc)) {
    tooltip_text <- paste0(
      tooltip_text,
      "<br>Description: ", df_volcano[[col_desc]]
    )
  }

  tooltip_text <- paste0(
    tooltip_text,
    "<br>log2FC: ", df_volcano$logFC,
    "<br>P.Value: ", df_volcano$P.Value
  )

  df_volcano$tooltip <- tooltip_text

  g_volcano <- ggplot(
    df_volcano,
    aes(
      x = logFC,
      y = -log10(P.Value),
      color = point_color,
      text = tooltip
    )
  ) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = colors_volcano) +
    labs(
      title = paste("Volcano plot:", comparison),
      x = "log2FC",
      y = "-log10(P.Value)"
    ) +
    theme(
      axis.ticks = element_line(),
      legend.position = "bottom",
      legend.key = element_rect(fill = "white"),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      text = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    )

  ############################ PROFILE PLOT #########################################

  expr <- as.data.frame(stats_data$data)
  annot <- as.data.frame(stats_data$annotation)
  meta  <- as.data.frame(stats_data$metadata)

  # Add phosphosite column
  expr$uniprot_id <- annot$uniprot_id

  # If no gene selected → return NULL for the profile plot
  if (is.null(user_phos_list) || user_phos_list == "") {
    return(list(
      volcano_plot = plotly::ggplotly(g_volcano, tooltip = "text"),
      profile_plot = NULL
    ))
  }

  # Filter to keep only the selected gene
  expr_sel <- expr[expr$uniprot_id == user_phos_list, , drop = FALSE]

  # Reshape to long
  df_profile <- expr_sel %>%
    tidyr::pivot_longer(cols = setdiff(colnames(expr_sel), "uniprot_id"),
                        names_to = "Sample_ID",
                        values_to = "LFQ_intensity") %>%
    dplyr::left_join(meta, by = c("Sample_ID" = "Simplify_ID"))

  df_profile$Sample_ID <- factor(
    df_profile$Sample_ID,
    levels = unique(df_profile$Sample_ID[order(df_profile$group)])
  )

  g_profile <- ggplot(
    df_profile,
    aes(
      x = Sample_ID,
      y = LFQ_intensity,
      group = uniprot_id,
      color = group,
      text = paste0("Phosphosite: ", uniprot_id)
    )
  ) +
    geom_line(size = 1.2) +
    geom_point(size = 2.2) +
    labs(
      title = paste("Intensity profile: ", user_phos_list),
      x = "Sample",
      y = "LFQ intensity (log2)"
    ) +
    theme(
      axis.ticks = element_line(),
      legend.position = "bottom",
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      text = element_text(size = 10),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10)
    )

  list(
    volcano_plot = plotly::ggplotly(g_volcano, tooltip = "text"),
    profile_plot = plotly::ggplotly(g_profile, tooltip = "text")
  )
}

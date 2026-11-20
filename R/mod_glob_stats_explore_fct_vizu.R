#' Generate volcano and protein profile plots
#'
#' @description
#' Generates visualizations to explore differential expression results.
#' This function produces:
#' \itemize{
#'   \item A volcano plot showing log fold-changes versus adjusted p-values
#'         with color-coded significance (up/down/non-significant) and optional highlighting
#'         of user-selected genes.
#'   \item A protein intensity profile plot for a selected gene, showing LFQ intensities
#'         across samples grouped by experimental condition.
#' }
#'
#' @param stats_data A list returned by the statistical analysis module, containing:
#'   \itemize{
#'     \item \code{results}: differential expression results per comparison
#'     \item \code{data}: quantitative data (LFQ intensities)
#'     \item \code{annotation}: protein annotation data
#'     \item \code{metadata}: experimental metadata for the samples
#'   }
#' @param comparison Character. Name of the comparison to plot (must match a key in \code{stats_data$results}).
#' @param user_gene_list Character vector of selected genes to highlight in the volcano plot and show in the profile plot. If NULL or empty, the profile plot is returned as NULL.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{volcano_plot}: a \pkg{plotly} object of the volcano plot
#'   \item \code{profile_plot}: a \pkg{plotly} object of the selected gene's intensity profile, or NULL if no gene is specified
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage
#' stats_data <- run_my_stats_module(my_data)
#' viz <- mod_glob_stats_explore_fct_vizu(stats_data, comparison = "CondA_vs_CondB", user_gene_list = "TP53")
#' viz$volcano_plot
#' viz$profile_plot
#' }
#'
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @importFrom dplyr left_join case_when
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @export

mod_glob_stats_explore_fct_vizu <- function(stats_data, comparison, user_gene_list = NULL) {

  if ("Protein.Group" %in% colnames(stats_data$annotation)) {
    col_pg <- "Protein.Group"
    col_gene <- "Genes"
    col_desc <- "First.Protein.Description"
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
  if (!is.null(user_gene_list) & user_gene_list != "") {
    df_volcano$point_color[df_volcano[[col_gene]] == user_gene_list] <- "selected"
  }

  colors_volcano <- c(
    "DownRegulated" = "#5555FF",
    "nonSignificant" = "grey",
    "UpRegulated" = "#FF5555",
    "selected" = "black"
  )

  g_volcano <- ggplot(df_volcano, aes(x = logFC, y = -log10(P.Value), color = point_color,
                                      text = paste0("Gene: ", .data[[col_gene]],
                                                    "<br>Description: ", .data[[col_desc]],
                                                    "<br>log2FC: ", logFC,
                                                    "<br>P.Value: ", P.Value))) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = colors_volcano) +
    labs(title = paste("Volcano plot: ", comparison),
         x = "log2FC",
         y = "-log10(P.Value)") +
    theme(axis.ticks = element_line(),
          legend.position = "bottom",
          legend.key = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10))

  ############################ PROFILE PLOT #########################################

  expr <- as.data.frame(stats_data$data)
  annot <- stats_data$annotation
  meta  <- stats_data$metadata

  # Add gene column
  expr[[col_gene]] <- annot[[col_gene]]

  # If no gene selected → return NULL for the profile plot
  if (is.null(user_gene_list) || user_gene_list == "") {
    return(list(
      volcano_plot = plotly::ggplotly(g_volcano, tooltip = "text"),
      profile_plot = NULL
    ))
  }

  # Filter to keep only the selected gene
  expr_sel <- expr[expr[[col_gene]] == user_gene_list, , drop = FALSE]

  # Reshape to long
  df_profile <- expr_sel %>%
    tidyr::pivot_longer(cols = setdiff(colnames(expr_sel), col_gene),
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
      group = .data[[col_gene]],
      color = group,
      text = paste0("Gene: ", .data[[col_gene]])
    )
  ) +
    geom_line(size = 1.2) +
    geom_point(size = 2.2) +
    labs(
      title = paste("Intensity profile: ", user_gene_list),
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

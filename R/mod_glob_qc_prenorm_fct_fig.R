#' Generate pre-normalization quality control plots
#'
#' @description
#' Generates a set of visualizations to assess the quality of quantitative proteomics
#' data before normalization. This function produces interactive \pkg{plotly} plots
#' and static visualizations that allow the user to examine intensity distributions,
#' variability, missing values, rank-abundance, correlations, and PCA structure.
#'
#' @details
#' The function expects a list returned by \code{mod_glob_qc_prenorm_utils_data_prep()}
#' (see \code{mod_glob_qc_prenorm_utils_data_prep.R}), which contains:
#' \itemize{
#'   \item reshaped intensity data in long format
#'   \item summary statistics per protein and per group
#'   \item counts of valid values (non-missing) per sample and per group
#'   \item rank-abundance and total intensity data
#'   \item correlation matrices and PCA results
#'   \item heatmap data
#' }
#'
#' The generated plots include:
#' \itemize{
#'   \item Density plots of log2 LFQ intensities per sample and per group
#'   \item CV (coefficient of variation) distribution per protein and per group
#'   \item Boxplots of log2 LFQ intensities per sample and per group
#'   \item Bar plots of the number of identified proteins and total input intensity
#'   \item Scatterplots of protein count vs total intensity
#'   \item Rank-abundance plots
#'   \item Sample correlation heatmaps
#'   \item PCA biplots
#'   \item Heatmaps of top variable proteins
#'   \item Distribution of proteins by number of valid values per sample and per group
#' }
#'
#' Most plots are returned as interactive \pkg{plotly} objects to allow tooltips
#' and dynamic exploration.
#'
#' @param qc_data A named list produced by \code{mod_glob_qc_prenorm_utils_data_prep()},
#'   containing preprocessed data and summary statistics required for plotting.
#'
#' @return A named list of plot objects, including interactive \pkg{plotly} plots,
#'   PCA biplots, and heatmaps. Elements can be rendered directly in Shiny or R Markdown.
#'
#' @examples
#' \dontrun{
#' # Prepare QC data
#' qc_data <- mod_glob_qc_prenorm_utils_data_prep(my_data)
#'
#' # Generate QC plots
#' qc_plots <- mod_glob_qc_prenorm_fct_fig(qc_data)
#'
#' # Display a density plot
#' qc_plots$s_distrib_int
#' }
#'
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @importFrom forcats fct_inorder
#' @importFrom scales scientific
#' @export

mod_glob_qc_prenorm_fct_fig <- function(qc_data) {

  if ("Protein.Group" %in% colnames(qc_data$data_long)) {
    col_int <- "Log2_LFQ_Intensity"
    col_pg <- "Protein.Group"
    col_gene <- "Genes"
  } else if ("PG.ProteinGroups" %in% colnames(qc_data$data_long)) {
    col_int <- "Shifted_Log2_LFQ_Intensity"
    col_pg <- "PG.ProteinGroups"
    col_gene <- "PG.Genes"
  }

  ############################ DENSITY PLOTS ##############################################

  ### SAMPLE INTENSITY (LOG2) DISTRIBUTION

  s_distrib_int <- ggplot(data = qc_data$data_long, mapping = aes(.data[[col_int]], color = Simplify_ID, text = Simplify_ID)) +
    stat_density(geom = 'line', position = "identity") +
    scale_x_continuous(name = "LFQ intensity (log2)") +
    scale_color_discrete(name = "Sample ID") +
    theme(axis.ticks = element_line(),
          legend.position = ifelse(length(unique(qc_data$data_long$Simplify_ID)) < 15, "bottom", "none"),
          legend.key = element_rect(fill="white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 10),
          legend.spacing.x = unit(0.25, 'cm'),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10))


  ### SAMPLE INTENSITY (LOG2) DISTRIBUTION BY GROUP
  g_distrib_int <- ggplot(data = qc_data$data_long, mapping = aes(.data[[col_int]], color = group, text = group)) +
    stat_density(geom = 'line', position = "identity") +
    scale_x_continuous(name = "LFQ intensity (log2)") +
    scale_color_discrete(name = "Condition") +
    theme(axis.ticks = element_line(),
          legend.position = ifelse(length(unique(qc_data$data_long$Simplify_ID)) < 15, "bottom", "none"),
          legend.key = element_rect(fill="white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 10),
          legend.spacing.x = unit(0.25, 'cm'),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10))

  ### CV DISTRIBUTION PER PROTEIN PER GROUP WITH 100 % VALID VALUES (AVOIR BIAIS IN MEAN CALCULATION)

  g_cv_mean_distrib <- ggplot() +
    geom_point(data = qc_data$stats_data, mapping = aes(x = Mean, y = CV, color = group,
                                                        text = paste("Protein Groups:", .data[[col_pg]],
                                                                     "<br>Genes:", .data[[col_gene]])),
               alpha = 0.5) +
    scale_x_continuous(name = "Mean Intensity") +
    scale_y_continuous(name = "CV (%)", limits = c(floor(min(qc_data$stats_data$CV)), ceiling(max(qc_data$stats_data$CV))),
                       breaks = c(seq(0,100,25), 200, 300, 400)) +
    theme(axis.ticks = element_line(),
          legend.position = "bottom",
          legend.key = element_rect(fill="white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 10),
          legend.spacing.x = unit(0.25, 'cm'),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10))

  g_cv_distrib <- ggplot(data = qc_data$stats_data, aes(x = group, y = CV, fill = group,
                                                        text = paste("Protein Groups:", .data[[col_pg]],
                                                                     "<br>Genes:", .data[[col_gene]]))) +
    geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA,
                 color = "black", size = 0.8, alpha = 0.7) +
    geom_point(data = filter(qc_data$stats_data, is_outlier),
               aes(x = group, y = CV,
                   text = paste("Protein Groups:", .data[[col_pg]],
                                "<br>Genes:", .data[[col_gene]])), size = 2) + # Add red points for outliers
    theme(axis.ticks = element_line(),
          legend.position = "bottom",
          legend.key = element_rect(fill="white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 10),
          legend.spacing.x = unit(0.25, 'cm'),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10))

  ############################## BOXPLOT ############################################
  s_bp_int  <- ggplot(data = qc_data$data_long, mapping = aes(x = forcats::fct_inorder(Simplify_ID), y = .data[[col_int]])) +
    geom_boxplot(aes(color = group),
                 position = position_dodge(0.9),
                 width = 0.2) +
    scale_y_continuous(name = "LFQ intensity (log2)") +
    theme(axis.ticks = element_line(),
          legend.position = "bottom",
          legend.key = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 10),
          legend.spacing.x = unit(0.25, 'cm'),
          axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 10))

  g_bp_int  <- ggplot(data = qc_data$data_long, mapping = aes(x = group, y = .data[[col_int]])) +
    geom_boxplot(aes(fill = group),
                 position = position_dodge(0.9),
                 width = 0.2) +
    scale_y_continuous(name = "LFQ intensity (log2)") +
    scale_x_discrete(name = "") +
    theme(axis.ticks = element_line(),
          legend.position = "bottom",
          legend.key = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 10),
          legend.spacing.x = unit(0.25, 'cm'),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10))

  #### Number of protein identified per sample

  s_id <- ggplot(data = qc_data$vv_data , mapping = aes(x = forcats::fct_inorder(Simplify_ID), y = sumVV, fill = group, text = sumVV)) +
    geom_bar(stat="identity") +
    scale_y_continuous(name = "Number of proteins identified") +
    scale_x_discrete(name = "Sample_ID") +
    theme(axis.ticks = element_line(),
          legend.position= "bottom",
          legend.key = element_rect(fill="white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 10),
          legend.spacing.x = unit(0.25, 'cm'),
          axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size=10))

  ### Bar plot with mean and sd per condition

  g_id  <- ggplot(data = qc_data$stats_vv_data) +
    geom_bar(aes(x = group, y = Mean, text = round(Mean)), stat="identity", fill="steelblue") +
    geom_errorbar(aes(x=group, ymin=Mean-SD, ymax=Mean+SD), width=0.25, colour="orange", linewidth=1.0) +
    scale_y_continuous(name = "Number of proteins identified (mean +/- sd)") +
    theme(axis.ticks = element_line(),
          legend.position= "bottom",
          legend.key = element_rect(fill="white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 10),
          legend.spacing.x = unit(0.25, 'cm'),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10))

  ### Bar plot of total input intensity per sample

  s_total_int <- ggplot(data = qc_data$sum_int_data) +
    geom_bar(mapping = aes(x = Simplify_ID, y = sumInt, fill=group, text = sumInt), stat = "identity") +
    scale_y_continuous(name = "Total input LFQ intensity") +
    theme(axis.ticks = element_line(),
          legend.position= "bottom",
          legend.key = element_rect(fill="white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 10),
          legend.spacing.x = unit(0.25, 'cm'),
          axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size=10))

  ###################### SCATTER PLOT ############################################

  #### Scatterplot between number of ID and Total input LFQ intensity

  id_vs_int <- ggplot() +
    geom_point(data = qc_data$ID_LFQ_data,
               mapping  = aes(x = sumVV, y = sumInt, colour = group, text = Simplify_ID),
               size = 3,
               shape=1) +
    scale_y_continuous(name = "Total input LFQ Intensity", labels = scales::scientific) +
    scale_x_continuous(name = "Total Number of Identified proteins") +
    theme(axis.ticks = element_line(),
          legend.position= "bottom",
          legend.title=element_blank(),
          legend.key = element_rect(fill="white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 14),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14))

  #### Rank Abundance per sample

  rank_plot <- ggplot() +
    geom_point(data = qc_data$rank_data,
               mapping  = aes(x = Rank, y = Mean,
                              text = paste("Protein Groups:", .data[[col_pg]],
                                           "<br>Genes:", .data[[col_gene]])),
               color = "black",
               size = 3,
               shape = 1,
               alpha = 0.5) +
    scale_y_continuous(name = "Mean LFQ intensity (log10)") +
    scale_x_continuous(name = "Protein rank abundance") +
    theme(axis.ticks = element_line(),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14))

  ########### CORRELATION PLOT ######################

  corr_plot <- ggcorrplot::ggcorrplot(corr = qc_data$corr_data[[1]],
                          method = "square",
                          hc.order = TRUE,
                          type = "lower",
                          outline.col = "white",
                          lab = TRUE,
                          p.mat = qc_data$corr_data[[2]],
                          lab_size = 2,
                          tl.cex = 8) +
    theme(plot.background = element_rect(fill = "white", color = NA))

  ################## PCA WITH PCAtools ##################################

  pca_plot <- PCAtools::biplot(qc_data$pca_data,
                     lab = qc_data$pca_data$metadata$Simplify_ID,
                     ntopLoadings = 5,
                     showLoadings = TRUE,
                     colby = 'group',
                     hline = 0, vline = 0,
                     legendPosition = 'right')

  ##### Heatmap of the top100 (for each PC1&PC2) most variable proteins from PCA

  col_side_colors_df <- data.frame(Condition = qc_data$meta_data$group)
  rownames(col_side_colors_df) <- rownames(qc_data$meta_data)

  heatmap_pcaload <- heatmaply::heatmaply(
    qc_data$heatmap_data,
    dist_method = "euclidean",
    hclust_method = "average",
    colors = c("blue", "white", "red"),
    col_side_colors = col_side_colors_df,
    key.title = "LFQ Intensity (log2)",
    plot_method = "plotly",  # ensures it returns a plotly object instead of saving PNG
    file = NULL              # prevents writing files
  )

  s_vv_distrib <- ggplot(qc_data$vv_distrib_sample, aes(x = valid_count, y = n, text = tooltip_text)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(
      title = "Protein Count by Number of Valid Values (All Samples)",
      x = "Number of Valid (non-NA) Values",
      y = "Number of Proteins"
    ) +
    theme(axis.ticks = element_line(),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14))

  g_vv_distrib <- ggplot(qc_data$vv_distrib_grp, aes(x = valid_count, y = n, fill = group, text = tooltip_text)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      title = "Protein Count by Number of Valid Values per Group",
      x = "Number of Valid (non-NA) Values",
      y = "Number of Proteins"
    ) +
    theme(axis.ticks = element_line(),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14))

  sample_annot <- data.frame(Condition = qc_data$meta_data$group)
  rownames(sample_annot) <- qc_data$meta_data$Simplify_ID

  if (nrow(qc_data$heatmap_on_off) > 0) {
    sample_annot <- data.frame(Condition = qc_data$meta_data$group)
    rownames(sample_annot) <- qc_data$meta_data$Simplify_ID

    heatmap_on_off <- heatmaply::heatmaply(
      qc_data$heatmap_on_off,
      dist_method = NULL,
      hclust_method = NULL,
      Rowv = FALSE,      # no row dendrogram
      Colv = FALSE,      # no column dendrogram
      colors = c("blue", "white", "red"),
      col_side_colors = sample_annot,
      key.title = "LFQ Intensity (log2)",
      plot_method = "plotly",  # ensures it returns a plotly object instead of saving PNG
      file = NULL              # prevents writing files
    )
  } else {
    heatmap_on_off <- NULL
  }


  list(
    s_distrib_int =  plotly::ggplotly(s_distrib_int, tooltip = "text"),
    g_distrib_int = plotly::ggplotly(g_distrib_int, tooltip = "text"),
    g_cv_mean_distrib = plotly::ggplotly(g_cv_mean_distrib, tooltip = "text"),
    g_cv_distrib = plotly::ggplotly(g_cv_distrib, tooltip = "text"),
    s_bp_int = plotly::ggplotly(s_bp_int, tooltip = "text"),
    g_bp_int = plotly::ggplotly(g_bp_int, tooltip = "text"),
    s_id = plotly::ggplotly(s_id, tooltip = "text"),
    g_id = plotly::ggplotly(g_id, tooltip = "text"),
    s_total_int = plotly::ggplotly(s_total_int, tooltip = "text"),
    id_vs_int = plotly::ggplotly(id_vs_int, tooltip = "text"),
    rank_plot = plotly::ggplotly(rank_plot, tooltip = "text"),
    corr_plot = plotly::ggplotly(corr_plot, tooltip = "text"),
    pca_plot = pca_plot,
    heatmap_pcaload = heatmap_pcaload,
    s_vv_distrib = plotly::ggplotly(s_vv_distrib, tooltip = "text"),
    g_vv_distrib = plotly::ggplotly(g_vv_distrib, tooltip = "text"),
    heatmap_on_off = heatmap_on_off
  )
}

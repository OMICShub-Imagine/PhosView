#' Generate pre-normalization QC plots for phosphosite data
#'
#' This function generates multiple QC visualizations for phosphosite intensity data
#' before normalization, including density plots, boxplots, CV distributions,
#' bar plots of identified phosphosites, total intensity plots, rank abundance plots,
#' correlation plots, PCA plots, and heatmaps.
#'
#' @param qc_data A list returned by \code{mod_glob_qc_prenorm_utils_data_prep()} (from
#'   \code{mod_glob_qc_prenorm_utils_data_prep.R}). This list should contain:
#'   \itemize{
#'     \item \code{data_long}: long-format phosphosite intensity data for density and boxplots
#'     \item \code{stats_data}: data frame with mean intensity, CV, and outlier flags
#'     \item \code{vv_data}: number of valid phosphosites per sample
#'     \item \code{stats_vv_data}: mean and SD per condition for valid phosphosites
#'     \item \code{sum_int_data}: total input intensity per sample
#'     \item \code{ID_LFQ_data}: scatterplot data of number of IDs vs. total LFQ intensity
#'     \item \code{rank_data}: protein rank abundance data
#'     \item \code{corr_data}: list of correlation matrices and p-values
#'     \item \code{pca_data}: PCAtools object with PCA results
#'     \item \code{heatmap_data}: matrix for PCA loadings heatmap
#'     \item \code{vv_distrib_sample}: distribution of valid values per sample
#'     \item \code{vv_distrib_grp}: distribution of valid values per group
#'     \item \code{meta_data}: metadata for samples (group, Simplify_ID, etc.)
#'     \item \code{heatmap_on_off}: optional on/off heatmap matrix
#'   }
#'
#' @return A named list of \pkg{ggplot2} and interactive \pkg{plotly} objects:
#' \describe{
#'   \item{s_distrib_int}{Density plot of LFQ intensities per sample}
#'   \item{g_distrib_int}{Density plot of LFQ intensities per group}
#'   \item{g_cv_mean_distrib}{Scatterplot of CV vs mean intensity per group}
#'   \item{g_cv_distrib}{Boxplot of CV per group with outliers highlighted}
#'   \item{s_bp_int}{Boxplot of LFQ intensities per sample}
#'   \item{g_bp_int}{Boxplot of LFQ intensities per group}
#'   \item{s_id}{Bar plot of phosphosites identified per sample}
#'   \item{g_id}{Bar plot of mean +/- SD of phosphosites identified per group}
#'   \item{s_total_int}{Bar plot of total LFQ intensity per sample}
#'   \item{id_vs_int}{Scatterplot of number of IDs vs total LFQ intensity per sample}
#'   \item{rank_plot}{Protein rank abundance plot}
#'   \item{corr_plot}{Correlation heatmap of samples}
#'   \item{pca_plot}{PCA biplot colored by group}
#'   \item{heatmap_pcaload}{Heatmap of top PCA loadings}
#'   \item{s_vv_distrib}{Distribution of valid phosphosites per sample}
#'   \item{g_vv_distrib}{Distribution of valid phosphosites per group}
#'   \item{heatmap_on_off}{Optional heatmap of on/off phosphosites (may be NULL)}
#' }
#'
#' @import ggplot2
#' @import plotly
#' @importFrom ggcorrplot ggcorrplot
#' @importFrom PCAtools biplot
#' @importFrom heatmaply heatmaply
#' @importFrom viridis viridis
#' @importFrom dplyr filter
#' @importFrom forcats fct_inorder
#' @export
#'
#' @examples
#' \dontrun{
#' qc_plots <- mod_phos_qc_prenorm_fct_fig(qc_data)
#' plotly::ggplotly(qc_plots$s_distrib_int)
#' }
mod_phos_qc_prenorm_fct_fig <- function(qc_data) {

  if ("Protein" %in% colnames(qc_data$data_long)) {
    col_int <- "Log2_LFQ_Intensity"
    col_pg <- "Protein"
    col_gene <- "Gene.Names"
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
                                                        text = paste("PTM.CollapseKey:", PTM.CollapseKey,
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

  g_cv_distrib <- ggplot(data = qc_data$stats_data, aes(x = group, y = CV, fill = group)) +
    geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA,
                 color = "black", size = 0.8, alpha = 0.7) +
    geom_point(data = filter(qc_data$stats_data, is_outlier),
               aes(x = group, y = CV,
                   text = paste("PTM.CollapseKey:", PTM.CollapseKey,
                                "<br>Genes:", .data[[col_gene]]),
                   size = 2)) + # Add red points for outliers
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
    scale_y_continuous(name = "Number of phosphosites identified") +
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
    scale_y_continuous(name = "Number of phosphosites identified (mean +/- sd)") +
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
    scale_x_continuous(name = "Total Number of Identified phosphosites") +
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
                              text = paste("PTM.CollapseKey:", PTM.CollapseKey,
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
      title = "Phosphosite Count by Number of Valid Values (All Samples)",
      x = "Number of Valid (non-NA) Values",
      y = "Number of Phosphosite"
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
      title = "Phosphosite Count by Number of Valid Values per Group",
      x = "Number of Valid (non-NA) Values",
      y = "Number of Phosphosite"
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

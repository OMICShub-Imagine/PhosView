#' Generate QC plots for phosphosite datasets
#'
#' This function generates a set of quality control (QC) plots for phosphosite
#' datasets at different filtering thresholds (unfiltered, 50%, 70%, 100% of replicates).
#' The plots summarize metrics such as mean LFQ intensity, localization probability,
#' coefficient of variation (CV), missing values, and replicate consistency.
#'
#' @param qc_data A \code{list} of DAList-like objects returned by
#'   \code{mod_phos_first_class_qc_utils_data_prep}. Expected elements:
#'   \describe{
#'     \item{phosphosite_unfilter}{unfiltered phosphosite dataset}
#'     \item{phosphosite_50pc}{dataset filtered at 50\% of replicates}
#'     \item{phosphosite_70pc}{dataset filtered at 70\% of replicates}
#'     \item{phosphosite_100pc}{dataset filtered at 100\% of replicates}
#'   }
#'
#' @details
#' The function produces multiple QC plots:
#' \itemize{
#'   \item \code{id}: barplot of mean valid values per sample with SD error bars
#'   \item \code{rank_plot}: phosphosite ranking based on mean LFQ intensity (log10) colored by mean localization probability
#'   \item \code{missing_vs_loc_boxplot}: boxplots of mean localization probability vs number of missing values
#'   \item \code{loc_vs_cv_plot}: scatterplot of localization probability vs CV per phosphosite
#'   \item \code{violin_cv_plot}: violin plot of CV per dataset
#'   \item \code{locprob_hist_plot}: histogram of mean localization probabilities
#'   \item \code{locprob_per_grp}: CV of localization probabilities per group
#'   \item \code{locprob_per_M}: CV of localization probabilities per group and multiplicity
#'   \item \code{locprob_LFQ_plot}: mean localization probability vs mean LFQ (log2)
#'   \item \code{CV_LFQ_plot}: CV vs mean LFQ intensity (log2) colored by mean localization probability
#' }
#'
#' @return A named \code{list} of \pkg{ggplot2} objects, one per QC plot.
#'
#' @import ggplot2
#' @importFrom stats cor na.omit quantile rnorm sd setNames
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import purrr
#' @import stringr
#' @export
#'
#' @examples
#' \dontrun{
#' # Prepare QC data
#' qc_data <- mod_phos_first_class_qc_utils_data_prep(
#'   phospho_unfilter = phospho_data,
#'   phospho_50pc = phospho_50pc,
#'   phospho_70pc = phospho_70pc,
#'   phospho_100pc = phospho_100pc,
#'   metadata = exp_metadata
#' )
#'
#' # Generate QC plots
#' qc_plots <- mod_phos_first_class_qc_fct_fig(qc_data)
#'
#' # View a plot
#' qc_plots$rank_plot
#' }

mod_phos_first_class_qc_fct_fig <- function(qc_data) {

  sample_metadata <- qc_data$phosphosite_unfilter$metadata

  phospho_site_unfilter <- dim(qc_data$phosphosite_unfilter$data)[1]

  site_id_phospho_unfilter <- data.frame(
    Sample = colnames(qc_data$phosphosite_unfilter$data),
    Valid_Values = colSums(!is.na(qc_data$phosphosite_unfilter$data))
  )

  phospho_site_100pc <- dim(qc_data$phosphosite_100pc$data)[1]

  site_id_phospho_100pc <- data.frame(
    Sample = colnames(qc_data$phosphosite_100pc$data),
    Valid_Values = colSums(!is.na(qc_data$phosphosite_100pc$data))
  )

  phospho_site_70pc <- dim(qc_data$phosphosite_70pc$data)[1]

  site_id_phospho_70pc <- data.frame(
    Sample = colnames(qc_data$phosphosite_70pc$data),
    Valid_Values = colSums(!is.na(qc_data$phosphosite_70pc$data))
  )

  phospho_site_50pc <- dim(qc_data$phosphosite_50pc$data)[1]

  site_id_phospho_50pc <- data.frame(
    Sample = colnames(qc_data$phosphosite_50pc$data),
    Valid_Values = colSums(!is.na(qc_data$phosphosite_50pc$data))
  )


  plot_data_sample_id <- rbind(site_id_phospho_unfilter,
                               site_id_phospho_50pc, site_id_phospho_70pc, site_id_phospho_100pc) %>%
    mutate(Category = c(rep("Phosphosite (unfiltered)", dim(sample_metadata)[1]),
                        rep("Phosphosite filtered at 50%", dim(sample_metadata)[1]),
                        rep("Phosphosite filtered at 70%", dim(sample_metadata)[1]),
                        rep("Phosphosite filtered at 100%", dim(sample_metadata)[1])))

  rownames(plot_data_sample_id) <- NULL

  plot_data_sample_id$Category <- factor(plot_data_sample_id$Category, levels = c("Phosphosite (unfiltered)",
                                                                                  "Phosphosite filtered at 50%",
                                                                                  "Phosphosite filtered at 70%",
                                                                                  "Phosphosite filtered at 100%"))

  summary_data <- plot_data_sample_id %>%
    group_by(Category) %>%
    summarise(
      mean_value = mean(Valid_Values, na.rm = TRUE),
      sd_value   = sd(Valid_Values, na.rm = TRUE),
      .groups = "drop"
    )

  # Plot with bars + SD error bars
  p1 <- ggplot(summary_data, aes(x = Category, y = mean_value, fill = Category)) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
                  width = 0.2,
                  position = position_dodge(width = 0.9)) +
    geom_text(aes(label = round(mean_value)),
              vjust = -1.5,
              position = position_dodge(width = 0.9),
              size = 3.5) +
    scale_fill_manual(values = c(
      "Phosphosite (unfiltered)"     = "#7CAE00",
      "Phosphosite filtered at 50%"     = "#00BFC4",
      "Phosphosite filtered at 70%" = "#C77CFF",
      "Phosphosite filtered at 100%" = "#F06060"
    )) +
    theme_minimal(base_size = 14) +
    scale_y_continuous(name = "Mean +/- SD of Valid Values",
                       limits = c(0, max(plot_data_sample_id$Valid_Values + 500)),
                       breaks = seq(0, max(plot_data_sample_id$Valid_Values + 500), 5000)) +
    scale_x_discrete(name = "") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  col_newname <- c(stringr::str_c("Locprob", sample_metadata$Simplify_ID, sep = "_"),
                   stringr::str_c("LFQ",  sample_metadata$Simplify_ID, sep = "_"))


  phosphosite_unfilter_quant <- qc_data$phosphosite_unfilter$data
  phosphosite_unfilter_locprob <- qc_data$phosphosite_unfilter$annotation %>%
    select(grep("PTM_SiteProbability", names(qc_data$phosphosite_unfilter$annotation), value = TRUE))
  phosphosite_unfilter_quant_prob <- cbind(phosphosite_unfilter_locprob, phosphosite_unfilter_quant)
  colnames(phosphosite_unfilter_quant_prob) <- col_newname

  phosphosite_50pc_quant <- qc_data$phosphosite_50pc$data
  phosphosite_50pc_locprob <- qc_data$phosphosite_50pc$annotation %>%
    select(grep("PTM_SiteProbability", names(qc_data$phosphosite_50pc$annotation), value = TRUE))
  phosphosite_50pc_quant_prob <- cbind(phosphosite_50pc_locprob, phosphosite_50pc_quant)
  colnames(phosphosite_50pc_quant_prob) <- col_newname

  phosphosite_70pc_quant <- qc_data$phosphosite_70pc$data
  phosphosite_70pc_locprob <- qc_data$phosphosite_70pc$annotation %>%
    select(grep("PTM_SiteProbability", names(qc_data$phosphosite_70pc$annotation), value = TRUE))
  phosphosite_70pc_quant_prob <- cbind(phosphosite_70pc_locprob, phosphosite_70pc_quant)
  colnames(phosphosite_70pc_quant_prob) <- col_newname

  phosphosite_100pc_quant <- qc_data$phosphosite_100pc$data
  phosphosite_100pc_locprob <- qc_data$phosphosite_100pc$annotation %>%
    select(grep("PTM_SiteProbability", names(qc_data$phosphosite_100pc$annotation), value = TRUE))
  phosphosite_100pc_quant_prob <- cbind(phosphosite_100pc_locprob, phosphosite_100pc_quant)
  colnames(phosphosite_100pc_quant_prob) <- col_newname

  # Function to process each dataframe
  process_df <- function(df){
    loccols <- grep('^Locprob_', colnames(df))
    lfqcols <- grep('^LFQ_', colnames(df))

    # compute metrics per row
    loc_mean <- rowMeans(df[, loccols, drop = FALSE], na.rm = TRUE)
    n_missing <- rowSums(is.na(df[, lfqcols, drop = FALSE]))
    lfq_mean <- rowMeans(df[, lfqcols, drop = FALSE], na.rm = TRUE)

    # 🔹 transform before taking mean
    lfq_log10 <- rowMeans(log10(df[, lfqcols, drop = FALSE]) + 10, na.rm = TRUE)
    lfq_log2  <- rowMeans(log2(df[, lfqcols, drop = FALSE]) + 10, na.rm = TRUE)

    lfq_sd <- apply(df[, lfqcols, drop = FALSE], 1, sd, na.rm = TRUE)
    cv <- ifelse(lfq_mean == 0, NA, (lfq_sd / lfq_mean) * 100)

    tibble(
      phosphosite = rownames(df),
      mean_locprob = loc_mean,
      n_missing = n_missing,
      mean_LFQ = lfq_mean,
      mean_LFQ_log10 = lfq_log10,
      mean_LFQ_log2 = lfq_log2,
      CV = cv
    )
  }


  # list of datasets
  datasets <- list(
    `Phosphosite (unfiltered)` = phosphosite_unfilter_quant_prob,
    `Phosphosite filtered at 50%` = phosphosite_50pc_quant_prob,
    `Phosphosite filtered at 70%` = phosphosite_70pc_quant_prob,
    `Phosphosite filtered at 100%` = phosphosite_100pc_quant_prob
  )

  processed <- purrr::imap_dfr(datasets, ~{
    df <- .x
    res <- process_df(df)
    res$dataset <- .y
    res
  }) %>% as.data.frame()


  processed$dataset <- factor(processed$dataset, levels = c("Phosphosite (unfiltered)",
                                                            "Phosphosite filtered at 50%",
                                                            "Phosphosite filtered at 70%",
                                                            "Phosphosite filtered at 100%"))

  # Ranking plot (X rows remove because NA in all sample in unfiltered dataset for X phosphositeS)
  # It means that they have locprob but are still considered as filtered or intensity <1....
  ranking_plot <- processed %>%
    group_by(dataset) %>%
    arrange(desc(mean_LFQ_log10)) %>%
    mutate(rank = row_number()+1) %>%
    ungroup() %>%
    ggplot(aes(x = rank, y = mean_LFQ_log10, color = mean_locprob)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_y_continuous(name = "Mean LFQ (log10)") +
    scale_x_continuous(name = "Protein rank") +
    scale_color_viridis_c(option = 'magma', name = 'Mean\nLocProb') +
    facet_wrap(~dataset, scales = 'free_x', ncol = 1) +
    labs(x = 'Rank (by mean LFQ)', y = 'log10(mean LFQ)') +
    theme_minimal()

  # Missing vs localisation probability
  # Compute counts per dataset x n_missing
  counts <- processed %>%
    group_by(dataset, n_missing) %>%
    summarise(n_sites = n(), .groups = "drop")

  # Violin/boxplot with counts annotated
  missing_vs_loc_boxplot <- processed %>%
    ggplot(aes(x = factor(n_missing), y = mean_locprob)) +
    geom_boxplot(fill = "skyblue", alpha = 0.6, outlier.size = 1.0, width = 0.6) +
    # add counts above each boxplot
    geom_text(
      data = counts,
      aes(
        x = factor(n_missing),
        y = 1.05,          # adjust y position above max
        label = n_sites
      ),
      inherit.aes = FALSE,
      size = 3
    ) +
    facet_wrap(~dataset, ncol = 2) +
    labs(
      x = "Number of missing LFQ values per phosphosite",
      y = "Mean localisation probability",
    ) +
    ylim(0, 1.1) +  # adjust to leave space for count annotations
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      strip.text = element_text(size = 12)  # facet labels
    )


  # Locprob vs CV
  loc_vs_cv_plot <- processed %>%
    ggplot(aes(x = mean_locprob, y = CV, color = mean_locprob)) +
    geom_point(alpha = 0.5, size = 0.8) +
    facet_wrap(~dataset, ncol = 2, scales = 'fixed') +
    scale_color_viridis_c(option = 'magma', name = 'Mean\nLocProb') +
    scale_y_continuous(name = "CV of raw LFQ intensity\nper phosphosite (%)", breaks = seq(0, 350, 25)) +
    scale_x_continuous(name = "Mean localisation probability") +
    theme_minimal()+
    theme(
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      strip.text = element_text(size = 12)  # facet labels
    )

  # Violin plot of CV
  # Removed 260 rows containing non-finite outside the scale range (because of non enough value to calculate CV)
  violin_cv_plot <- processed %>%
    ggplot(aes(x = dataset, y = CV)) +
    geom_violin(trim = TRUE, alpha = 0.6, fill = "skyblue") +
    geom_boxplot(width = 0.08, outlier.size = 1.0, alpha = 0.8) +
    scale_y_continuous(name = "CV of raw LFQ intensity\nper phosphosite (%)",
                       breaks = seq(0, 350, 25)) +
    scale_x_discrete(name = "") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      strip.text = element_text(size = 12)  # facet labels
    )

  locprob_hist_plot <- processed %>%
    ggplot(aes(x = mean_locprob, fill = dataset)) +
    geom_histogram(
      bins = 30,
      color = "black",   # black edges
      alpha = 0.6,
      position = "identity"
    ) +
    facet_wrap(~dataset, ncol = 2, scales = "fixed") +
    labs(
      x = "Mean localisation probability",
      y = "Count",
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.position = "none"
    )



  #### Phosphosite localisation probability distribution per group and per filtering method
  colnames(phosphosite_unfilter_locprob) <- col_newname[1:dim(sample_metadata)[1]]
  colnames(phosphosite_50pc_locprob) <- col_newname[1:dim(sample_metadata)[1]]
  colnames(phosphosite_70pc_locprob) <- col_newname[1:dim(sample_metadata)[1]]
  colnames(phosphosite_100pc_locprob) <- col_newname[1:dim(sample_metadata)[1]]

  # Function to compute CV per row for each group
  # Generic CV computation using metadata
  compute_cv <- function(df, dataset_name, sample_metadata, prefix = "Locprob_") {

    # Extract sample column names (those present in metadata)
    sample_cols <- names(df)[names(df) %in% paste0(prefix, sample_metadata$Simplify_ID)]

    # Match them to metadata to know which group each belongs to
    meta <- sample_metadata %>%
      filter(paste0(prefix, Simplify_ID) %in% names(df)) %>%
      mutate(Column = paste0(prefix, Simplify_ID))

    # Compute CV per group based on metadata grouping
    cv_data <- meta %>%
      split(.$group) %>%
      purrr::map_dfr(function(gmeta) {
        cols <- gmeta$Column
        cv <- apply(df[, cols, drop = FALSE], 1, function(x) (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100)
        tibble(
          Phosphosite = rownames(df),
          Group = unique(gmeta$group),
          CV = cv
        )
      })

    # Add metadata and derived columns
    cv_data %>%
      mutate(
        Dataset = dataset_name,
        Multiplicity = stringr::str_extract(Phosphosite, "_M[0-9]+"),
        Multiplicity = stringr::str_remove(Multiplicity, "_"),
        Group_Multiplicity = stringr::str_c(Group, Multiplicity, sep = "_")
      )
  }

  cv_data <- bind_rows(
    compute_cv(phosphosite_unfilter_locprob, "Phosphosite (unfiltered)", sample_metadata),
    compute_cv(phosphosite_50pc_locprob, "Phosphosite filtered at 50%", sample_metadata),
    compute_cv(phosphosite_70pc_locprob, "Phosphosite filtered at 70%", sample_metadata),
    compute_cv(phosphosite_100pc_locprob, "Phosphosite filtered at 100%", sample_metadata)
  )

  cv_data$Dataset <- factor(
    cv_data$Dataset,
    levels = c(
      "Phosphosite (unfiltered)",
      "Phosphosite filtered at 50%",
      "Phosphosite filtered at 70%",
      "Phosphosite filtered at 100%"
    )
  )

  # Plot CV distribution per dataset and group
  locprob_per_grp <- ggplot(cv_data, aes(x = Group, y = CV, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.4) +
    theme_minimal(base_size = 14) +
    scale_y_continuous(name = "CV of localisation probability (%)", breaks = seq(0, 300, 25)) +
    scale_x_discrete(name = "", breaks = NULL) +
    facet_wrap(~ Dataset, scales = "fixed") +
    scale_fill_brewer(palette = "Set1")

  locprob_per_M <- ggplot(cv_data, aes(x = Group_Multiplicity, y = CV, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.4) +
    theme_minimal(base_size = 14) +
    scale_y_continuous(name = "CV of localisation probability (%)", breaks = seq(0, 300, 25)) +
    scale_x_discrete(name = "", labels = c("CTRL_M1" = "M1",
                                           "CTRL_M2" = "M2",
                                           "CTRL_M3" = "M3",
                                           "EGF_M1" = "M1",
                                           "EGF_M2" = "M2",
                                           "EGF_M3" = "M3")) +
    facet_wrap(~ Dataset, scales = "fixed") +
    scale_fill_brewer(palette = "Set1")


  #### Mean localisation probability against mean log2 LFQ intensities per phosphosite

  locprob_LFQ_plot <- processed %>%
    ggplot(aes(x = mean_locprob, y = mean_LFQ_log2, color = CV)) +
    geom_point(
      alpha = 0.6,
      position = "identity"
    ) +
    scale_color_viridis_c(
      option = "magma",
      name = "CV",
      limits = c(0, 100),
      oob = scales::squish
    ) +
    facet_wrap(~dataset, ncol = 2, scales = "fixed") +
    labs(
      x = "Mean localisation probability",
      y = "Mean LFQ intensity (log2)",
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      strip.text = element_text(size = 12),
    )



  CV_LFQ_plot <- processed %>%
    ggplot(aes(x = mean_LFQ_log2, y = CV, color = mean_locprob)) +
    geom_point(
      alpha = 0.6,
      position = "identity"
    ) +
    scale_color_viridis_c(
      option = "magma",
      name = "mean_locprob",
    ) +
    facet_wrap(~dataset, ncol = 2, scales = "fixed") +
    scale_y_continuous(name = "CV" , breaks = seq(0, 300, 50)) +
    scale_x_continuous(name = "Mean LFQ intensity (log2)", breaks = seq(10, 25, 2.5)) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      strip.text = element_text(size = 12),
    )





  list(
    id = p1,
    rank_plot = ranking_plot,
    missing_vs_loc_boxplot = missing_vs_loc_boxplot,
    loc_vs_cv_plot = loc_vs_cv_plot,
    violin_cv_plot = violin_cv_plot,
    locprob_hist_plot = locprob_hist_plot,
    locprob_per_grp = locprob_per_grp,
    locprob_per_M = locprob_per_M,
    locprob_LFQ_plot = locprob_LFQ_plot,
    CV_LFQ_plot = CV_LFQ_plot
  )

}

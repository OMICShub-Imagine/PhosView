#' Volcano plot for Global Analysis
#'
#' @description
#' Generates an interactive volcano plot for global proteomics analysis,
#' displaying log fold-change versus adjusted p-values for a selected comparison.
#' Points are color-coded according to differential expression status
#' (up-regulated, down-regulated, non-significant), with optional highlighting
#' of user-selected genes.
#'
#' @param data A data.frame or tibble containing processed global analysis results.
#'   Must include columns corresponding to the selected comparison, such as
#'   log fold-change, adjusted p-values, and significance flags.
#' @param comparison Character scalar specifying the comparison name used to
#'   construct column names (e.g. \code{"ConditionA_vs_ConditionB"}).
#' @param user_gene_list Optional character vector of gene names to highlight
#'   in the volcano plot.
#'
#' @return A \pkg{plotly} object representing the interactive volcano plot.
#'
#' @importFrom plotly plot_ly layout
#' @export
renderVolcano_Glob <- function(data, comparison, user_gene_list = NULL) {
  data$point_color <- ifelse(
    data[[paste0(comparison, ".sig.FDR_GlobalAnalysis")]] == -1, "DownRegulated",
    ifelse(data[[paste0(comparison, ".sig.FDR_GlobalAnalysis")]] == 1, "UpRegulated", "nonSignificant")
  )

  colors <- c("DownRegulated" = "#5555FF", "nonSignificant" = "grey", "UpRegulated" = "#FF5555")
  if (!is.null(user_gene_list) && length(user_gene_list) > 0) {
    data$point_color[data$PG.Genes %in% user_gene_list] <- "selected"
    colors <- c(colors, "selected" = "black")
  }

  plotly::plot_ly(
    data,
    x = ~get(paste0(comparison, ".logFC_GlobalAnalysis")),
    y = ~-log10(get(paste0(comparison, ".adj.P.Val_GlobalAnalysis"))),
    type = 'scatter',
    mode = 'markers',
    color = ~point_color,
    colors = colors,
    marker = list(opacity = 0.5),
    text = ~paste0(
      "Gene: ", get("PG.Genes"), "<br>",
      "Description: ", get("PG.ProteinDescriptions"), "<br>",
      "LogFC: ", get(paste0(comparison, ".logFC_GlobalAnalysis")), "<br>",
      "-log10(p-value): ", -log10(get(paste0(comparison, ".adj.P.Val_GlobalAnalysis")))
    ),
    hoverinfo = 'text'
  ) %>%
    plotly::layout(
      title = "Volcano Plot (Global Analysis)",
      xaxis = list(title = paste0(comparison, ".logFC")),
      yaxis = list(title = paste0("-log10(", comparison, ".adj.P.Val)"))
    )
}

#' Volcano plot for Phospho Analysis
#'
#' @description
#' Generates an interactive volcano plot for phosphoproteomics analysis.
#' The plot displays log fold-change versus adjusted p-values for a selected
#' comparison, with points colored according to phosphosite-level differential
#' regulation and optional highlighting of selected genes.
#'
#' @param data A data.frame or tibble containing processed phosphoproteomics
#'   analysis results, including phosphosite annotation.
#' @param comparison Character scalar specifying the comparison name used to
#'   construct column names.
#' @param user_gene_list Optional character vector of gene names to highlight
#'   in the volcano plot.
#'
#' @return A \pkg{plotly} object representing the interactive phospho volcano plot.
#'
#' @importFrom plotly plot_ly layout
#' @export
renderVolcano_Phos <- function(data, comparison, user_gene_list = NULL) {
  data$point_color <- ifelse(
    data[[paste0(comparison, ".sig.FDR_PhosphoAnalysis")]] == -1, "DownRegulated",
    ifelse(data[[paste0(comparison, ".sig.FDR_PhosphoAnalysis")]] == 1, "UpRegulated", "nonSignificant")
  )

  colors <- c("DownRegulated" = "#5555FF", "nonSignificant" = "grey", "UpRegulated" = "#FF5555")
  if (!is.null(user_gene_list) && length(user_gene_list) > 0) {
    data$point_color[data$PG.Genes %in% user_gene_list] <- "selected"
    colors <- c(colors, "selected" = "black")
  }

  plotly::plot_ly(
    data,
    x = ~get(paste0(comparison, ".logFC_PhosphoAnalysis")),
    y = ~-log10(get(paste0(comparison, ".adj.P.Val_PhosphoAnalysis"))),
    type = 'scatter',
    mode = 'markers',
    color = ~point_color,
    colors = colors,
    marker = list(opacity = 0.5),
    text = ~paste0(
      "PhosphoSite: ", get("PTM.CollapseKey"), "<br>",
      "Gene: ", get("PG.Genes"), "<br>",
      "Description: ", get("PG.ProteinDescriptions"), "<br>",
      "LogFC: ", get(paste0(comparison, ".logFC_PhosphoAnalysis")), "<br>",
      "-log10(p-value): ", -log10(get(paste0(comparison, ".adj.P.Val_PhosphoAnalysis")))
    ),
    hoverinfo = 'text'
  ) %>%
    plotly::layout(
      title = "Volcano Plot (Phospho Analysis)",
      xaxis = list(title = paste0(comparison, ".logFC")),
      yaxis = list(title = paste0("-log10(", comparison, ".adj.P.Val)"))
    )
}

#' Protein intensity profile plot (Global Analysis)
#'
#' @description
#' Generates an interactive protein intensity profile plot for global analysis.
#' The plot shows LFQ intensities across samples, ordered by experimental group,
#' for selected genes.
#'
#' @param data A data.frame or tibble containing global LFQ intensity values
#'   and protein annotations.
#' @param metadata A data.frame containing sample metadata, including sample
#'   identifiers and experimental group assignments.
#' @param user_gene_list Optional character vector of gene names to display
#'   in the profile plot.
#'
#' @return A \pkg{plotly} object representing the global protein intensity profile.
#'
#' @importFrom plotly plot_ly layout
#' @importFrom dplyr select mutate
#' @importFrom tidyr pivot_longer
#' @export
renderPP_Glob <- function(
    data,
    metadata,
    user_gene_list = NULL,
    gene_col_glob,
    protein_id_col_glob
) {

  gene_col <- paste0(gene_col_glob, "_GlobalAnalysis")
  protein_col <- paste0(protein_id_col_glob, "_GlobalAnalysis")

  df_pp <- data %>%
    dplyr::select(dplyr::any_of(c(
      paste(metadata$Simplify_ID, "GlobalAnalysis", sep = "_"),
      gene_col,
      protein_col
    ))) %>%
    tidyr::pivot_longer(
      cols = 1:length(metadata$Simplify_ID),
      names_to = "Sample_ID",
      values_to = "LFQ_intensity"
    ) %>%
    dplyr::mutate(Group = rep(metadata$group, nrow(data))) %>%
    dplyr::mutate(Sample_ID = stringr::str_remove(Sample_ID, '_GlobalAnalysis'))

  df_pp$Sample_ID <- factor(df_pp$Sample_ID,
                            levels = unique(df_pp$Sample_ID[order(df_pp$Group)]))

  plotly::plot_ly(
    subset(df_pp, df_pp[[gene_col]] %in% user_gene_list),
    x = ~Sample_ID,
    y = ~LFQ_intensity,
    type = 'scatter',
    mode = 'lines+markers',
    color = ~forcats::fct_inorder(Group),
    text = ~paste0("Gene: ", get(gene_col)),
    hoverinfo = 'text'
  )
}

#' Phosphosite intensity profile plot (Phospho Analysis)
#'
#' @description
#' Generates an interactive phosphosite intensity profile plot showing LFQ
#' intensities across samples for selected genes. Multiple phosphosites
#' belonging to the same protein can be visualized simultaneously.
#'
#' @param data A data.frame or tibble containing phosphoproteomics LFQ intensity
#'   values and phosphosite annotations.
#' @param metadata A data.frame containing sample metadata, including sample
#'   identifiers and experimental group assignments.
#' @param user_gene_list Optional character vector of gene names to display
#'   in the profile plot.
#'
#' @return A \pkg{plotly} object representing the phosphosite intensity profile.
#'
#' @importFrom plotly plot_ly layout
#' @importFrom dplyr select mutate
#' @importFrom tidyr pivot_longer
#' @export
renderPP_Phos <- function(
    data,
    metadata,
    user_gene_list = NULL,
    gene_col_phos,
    protein_id_col_phos,
    phospho_site_col
) {

  gene_col <- paste0(gene_col_phos, "_PhosphoAnalysis")
  protein_col <- paste0(protein_id_col_phos, "_PhosphoAnalysis")
  site_col <- paste0(phospho_site_col, "_PhosphoAnalysis")

  df_pp <- data %>%
    dplyr::select(dplyr::any_of(c(
      paste(metadata$Simplify_ID, "PhosphoAnalysis", sep = "_"),
      gene_col,
      protein_col,
      site_col
    ))) %>%
    tidyr::pivot_longer(
      cols = 1:length(metadata$Simplify_ID),
      names_to = "Sample_ID",
      values_to = "LFQ_intensity"
    ) %>%
    dplyr::mutate(Group = rep(metadata$group, nrow(data))) %>%
    dplyr::mutate(Sample_ID = stringr::str_remove(Sample_ID, '_PhosphoAnalysis'))

  df_pp$Sample_ID <- factor(df_pp$Sample_ID,
                            levels = unique(df_pp$Sample_ID[order(df_pp$Group)]))

  plotly::plot_ly(
    subset(df_pp, df_pp[[gene_col]] %in% user_gene_list),
    x = ~Sample_ID,
    y = ~LFQ_intensity,
    type = 'scatter',
    mode = 'lines+markers',
    color = ~forcats::fct_inorder(Group),
    split = ~get(site_col),
    text = ~paste0(
      "Gene: ", get(gene_col),
      "<br>PhosphoSite: ", get(site_col)
    ),
    hoverinfo = 'text'
  )
}

#' Dual fold-change (Global vs Phospho) plot
#'
#' @description
#' Generates an interactive dual fold-change scatter plot comparing global
#' protein-level and phosphosite-level log fold-changes for a given comparison.
#' Points are color-coded according to combined differential regulation patterns
#' between global and phospho analyses, with optional highlighting of selected genes.
#'
#' @param global_data A data.frame or tibble containing global analysis results,
#'   including log fold-changes and significance flags.
#' @param phospho_data A data.frame or tibble containing phosphoproteomics analysis
#'   results, including log fold-changes and significance flags.
#' @param comparison Character scalar specifying the comparison name used to
#'   construct column names.
#' @param user_gene_list Optional character vector of gene names to highlight
#'   in the dual fold-change plot.
#'
#' @return A \pkg{plotly} object representing the dual fold-change scatter plot.
#'
#' @importFrom plotly plot_ly layout
#' @importFrom dplyr select mutate case_when
#' @export
renderFC <- function(
    global_data,
    phospho_data,
    comparison,
    protein_id_col_glob,
    protein_id_col_phos,
    gene_col_glob,
    gene_col_phos,
    desc_col_phos,
    phospho_site_col,
    user_gene_list = NULL
) {

  logfc_glob_col <- paste0(comparison, ".logFC_GlobalAnalysis")
  sig_glob_col   <- paste0(comparison, ".sig.FDR_GlobalAnalysis")

  logfc_phos_col <- paste0(comparison, ".logFC_PhosphoAnalysis")
  sig_phos_col   <- paste0(comparison, ".sig.FDR_PhosphoAnalysis")

  protein_glob_col <- paste0(protein_id_col_glob, "_GlobalAnalysis")
  protein_phos_col <- paste0(protein_id_col_phos, "_PhosphoAnalysis")

  gene_glob_col <- paste0(gene_col_glob, "_GlobalAnalysis")
  gene_phos_col <- paste0(gene_col_phos, "_PhosphoAnalysis")

  desc_phos_col <- if(!is.null(desc_col_phos)) paste0(desc_col_phos, "_PhosphoAnalysis") else NULL
  site_col <- paste0(phospho_site_col, "_PhosphoAnalysis")

  FC_glob <- global_data %>%
    dplyr::select(dplyr::any_of(c(
      logfc_glob_col,
      sig_glob_col,
      protein_glob_col,
      gene_glob_col
    )))

  FC_phos <- phospho_data %>%
    dplyr::select(dplyr::any_of(c(
      logfc_phos_col,
      sig_phos_col,
      protein_phos_col,
      gene_phos_col,
      desc_phos_col,
      site_col
    )))

  FC_fused <- base::merge(
    FC_glob,
    FC_phos,
    by.x = protein_glob_col,
    by.y = protein_phos_col,
    all = FALSE
  )

  FC_fused$Identifier <- dplyr::coalesce(
    FC_fused[[gene_glob_col]],
    FC_fused[[gene_phos_col]]
  )

  FC_fused <- FC_fused %>%
    dplyr::mutate(point_color = dplyr::case_when(
      (.[[sig_phos_col]] == -1 & .[[sig_glob_col]] == -1) ~ "Both DownReg",
      (.[[sig_phos_col]] == 1 & .[[sig_glob_col]] == 1) ~ "Both UpReg",
      (.[[sig_phos_col]] == -1 & .[[sig_glob_col]] == 1) ~ "DownReg in Phospho & UpReg in Global",
      (.[[sig_phos_col]] == 1 & .[[sig_glob_col]] == -1) ~ "UpReg in Phospho & DownReg in Global",
      (.[[sig_phos_col]] == 1 & .[[sig_glob_col]] == 0) ~ "Upreg in Phospho & NonSign in Global",
      (.[[sig_phos_col]] == -1 & .[[sig_glob_col]] == 0) ~ "Downreg in Phospho & NonSign in Global",
      (.[[sig_phos_col]] == 0 & .[[sig_glob_col]] == 1) ~ "NonSign in Phospho and Upreg in Global",
      (.[[sig_phos_col]] == 0 & .[[sig_glob_col]] == -1) ~ "NonSign in Phospho and Downreg in Global",
      TRUE ~ "NonSign in both"
    ))

  colors <- c(
    "Both DownReg" = "#5555FF",
    "Both UpReg" = "#FF5555",
    "DownReg in Phospho & UpReg in Global" = "#00AA00",
    "UpReg in Phospho & DownReg in Global" = "#AA00AA",
    "Upreg in Phospho & NonSign in Global" = "#FFAA00",
    "Downreg in Phospho & NonSign in Global" = "#55FF55",
    "NonSign in Phospho and Upreg in Global" = "#0000AA",
    "NonSign in Phospho and Downreg in Global" = "#AA0000",
    "NonSign in both" = "grey"
  )

  if (!is.null(user_gene_list) && length(user_gene_list) > 0) {
    FC_fused$point_color[FC_fused$Identifier %in% user_gene_list] <- "selected"
    colors <- c(colors, "selected" = "black")
  }

  plotly::plot_ly(
    FC_fused,
    x = ~get(logfc_glob_col),
    y = ~get(logfc_phos_col),
    type = 'scatter',
    mode = 'markers',
    color = ~point_color,
    colors = colors,
    marker = list(opacity = 0.6),
    text = ~paste0(
      "PhosphoSite: ", get(site_col), "<br>",
      "Gene: ", Identifier, "<br>",
      "Description: ", if(!is.null(desc_phos_col)) get(desc_phos_col) else NA, "<br>",
      "LogFC (Global): ", get(logfc_glob_col), "<br>",
      "LogFC (Phospho): ", get(logfc_phos_col), "<br>",
      "Significance (Global): ", get(sig_glob_col), "<br>",
      "Significance (Phospho): ", get(sig_phos_col)
    ),
    hoverinfo = 'text'
  ) %>%
    plotly::layout(
      title = "Dual-FC Plot",
      xaxis = list(title = logfc_glob_col),
      yaxis = list(title = logfc_phos_col)
    )
}

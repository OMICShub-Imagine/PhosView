#' Prepare QC data for plotting
#' @param quant_data should corresponds to data of \code{\link[proteoDA]{DAList}}
#' @param meta_data  should corresponds to metadata of \code{\link[proteoDA]{DAList}}
#' @param prot_annot should corresponds to annotation of \code{\link[proteoDA]{DAList}}
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom stats cor na.omit quantile rnorm sd setNames
#' @return A list of prepared data for QC plots

mod_glob_qc_postnorm_utils_data_prep <- function(quant_data, meta_data, prot_annot) {
  # Ensure Sample_ID columns exist in quant_data
  sample_cols <- intersect(meta_data$Simplify_ID, colnames(quant_data))
  if(length(sample_cols) == 0) stop("No matching sample columns found in quant data")

  if ("Protein.Group" %in% colnames(prot_annot)) {
    col_pg <- "Protein.Group"
    col_gene <- "Genes"
  } else if ("PG.ProteinGroups" %in% colnames(prot_annot)) {
    col_pg <- "PG.ProteinGroups"
    col_gene <- "PG.Genes"
  }

  # For plots that use data in long format
  data_long <- quant_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column(col_pg) %>%
    pivot_longer(cols = 2:(dim(quant_data)[2]+1),
                 names_to = "Simplify_ID",
                 values_to = "Log2_LFQ_Intensity") %>%
    left_join(prot_annot, by = col_pg) %>%
    left_join(meta_data, by = "Simplify_ID") %>%
    dplyr::mutate(group = as.factor(.$group)) %>%
    dplyr::mutate(Simplify_ID = factor(.$Simplify_ID, levels = unique(.$Simplify_ID[order(.$group)])))

  # Use for number of valid values per sample and per groups
  vv <- data_long %>%
    group_by(Simplify_ID, group) %>%
    summarise(sumVV = sum(!is.na(Log2_LFQ_Intensity))) %>%
    ungroup() %>%
    as.data.frame() %>%
    dplyr::mutate(Sample_ID = factor(Simplify_ID, levels = unique(Simplify_ID[order(group)])))

  stats_vv <- vv %>%
    group_by(group) %>%
    summarise(Mean = mean(sumVV), SD = sd(sumVV))


  # Use for the ranking plot
  rank <- quant_data %>%
    as.data.frame() %>%
    dplyr::mutate(VV_Number = rowSums(!is.na(.))) %>%
    filter(VV_Number >= 1) %>%
    dplyr::select(-VV_Number) %>%
    as.data.frame() %>%
    dplyr::mutate(Mean = rowMeans(., na.rm = TRUE)) %>%
    arrange(., desc(Mean)) %>%
    dplyr::mutate(Rank = 1:nrow(.)) %>%
    dplyr::select(c(Mean, Rank)) %>%
    rownames_to_column(col_pg) %>%
    left_join(prot_annot, by = col_pg)

  # Compute a correlation matrix
  corr <- round(cor(quant_data %>% as.data.frame() %>% na.omit()), 2)

  # Compute a matrix of correlation p-values
  p.mat <- ggcorrplot::cor_pmat(quant_data %>% as.data.frame() %>% na.omit())

  # Precomputation of PCA metrics
  pca_df <- quant_data %>%
    na.omit() %>%
    as.data.frame() %>%
    rownames_to_column(col_pg) %>%
    left_join(prot_annot, by = col_pg) %>%
    dplyr::mutate(Unique_ID = paste("PG:", .data[[col_pg]], "Gene:", .data[[col_gene]])) %>%
    column_to_rownames("Unique_ID") %>%
    select(!colnames(prot_annot))

  p <- PCAtools::pca(pca_df,
                     metadata = meta_data,
                     center = TRUE,
                     scale = TRUE)


  ### Input for heatmap

  # 1. Get unique groups, sorted alphabetically for consistency
  unique_groups <- sort(unique(meta_data$group))

  # 2. Define a fixed color palette (max 12 colors)
  max_colors <- 12
  palette_colors <- RColorBrewer::brewer.pal(n = max_colors, name = "Set3")

  # 3. Assign colors to groups in a reproducible way
  group_colors <- setNames(palette_colors[1:length(unique_groups)], unique_groups)

  # 4. Create annotation data frame for sample conditions
  annotation_df <- data.frame(Condition = meta_data$group)
  rownames(annotation_df) <- colnames(quant_data)  # Make sure row names match sample names

  # 5. Extract absolute values of loadings for PC1 and PC2
  abs_loadings <- abs(p$loadings[, c("PC1", "PC2")])

  # 6. Rank and extract the top 100 for each PC
  top_PC1 <- rownames(abs_loadings)[order(abs_loadings[, "PC1"], decreasing = TRUE)][1:100]
  top_PC2 <- rownames(abs_loadings)[order(abs_loadings[, "PC2"], decreasing = TRUE)][1:100]

  # 7. Combine and get unique gene names
  top_genes <- unique(c(top_PC1, top_PC2))

  # 8.Create a unique ID by combining PG.Genes and PG.ProteinGroups & Process the matrix (log2 transform, remove NAs scale)
  processed_matrix <- quant_data %>%
    t() %>%
    scale(center = TRUE, scale = TRUE) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(col_pg) %>%
    left_join(prot_annot, by = col_pg) %>%
    dplyr::mutate(Unique_ID = paste("PG:", .data[[col_pg]], "Gene:", .data[[col_gene]])) %>%
    column_to_rownames("Unique_ID") %>%
    select(!colnames(prot_annot)) %>%
    filter(row.names(.) %in% top_genes)

  # Some predefined dataframe for some plot about valid value distribution per sample

  # Total number of samples (used for range)
  n_samples <- data_long %>%
    pull(Simplify_ID) %>%
    n_distinct()

  # Count valid values per protein
  valid_counts_all <- data_long %>%
    group_by(uniprot_id) %>%
    summarise(valid_count = sum(!is.na(Log2_LFQ_Intensity)), .groups = "drop") %>%
    count(valid_count)

  # Fill in all values 0 to n_samples
  valid_counts_all_complete <- tibble(valid_count = 1:n_samples) %>%
    left_join(valid_counts_all, by = "valid_count") %>%
    dplyr::mutate(n = replace_na(n, 0))

  # Add % column
  total_proteins <- sum(valid_counts_all_complete$n)
  valid_counts_all_complete <- valid_counts_all_complete %>%
    dplyr::mutate(
      percent = round((n / total_proteins) * 100, 1),
      valid_count = factor(valid_count, levels = as.character(0:n_samples)),
      tooltip_text = paste0(
        "Valid Values: ", valid_count,
        "<br>Proteins: ", n,
        "<br>% of Total: ", percent, "%"
      )
    )

  # Some predefined dataframe for some plot about valid value distribution per group

  # Get total number of samples per group
  samples_per_group <- data_long %>%
    distinct(Simplify_ID, group) %>%
    count(group, name = "n_samples")

  # Count valid values per protein per group
  valid_counts_group <- data_long %>%
    group_by(group, uniprot_id) %>%
    summarise(valid_count = sum(!is.na(Log2_LFQ_Intensity)), .groups = "drop")

  # Tabulate counts of valid values per group
  distribution_group <- valid_counts_group %>%
    count(group, valid_count)

  # Get max number of samples across all groups
  max_samples <- samples_per_group$n_samples %>% max()

  # Create full grid of (group x 0:max_samples) to ensure all bars appear
  full_grid <- expand.grid(
    group = levels(data_long$group),
    valid_count = 0:max_samples
  )

  # Merge and fill missing with 0
  distribution_group_complete <- full_grid %>%
    left_join(distribution_group, by = c("group", "valid_count")) %>%
    dplyr::mutate(n = replace_na(n, 0))

  # Add total proteins per group for percentage calculation
  distribution_group_complete <- distribution_group_complete %>%
    group_by(group) %>%
    dplyr::mutate(
      total_proteins = sum(n),
      percent = round(100 * n / total_proteins, 1),
      valid_count = factor(valid_count, levels = as.character(0:max_samples)),
      tooltip_text = paste0(
        "Group: ", group,
        "<br>Valid Values: ", valid_count,
        "<br>Proteins: ", n,
        "<br>% of Group: ", percent, "%"
      )
    ) %>%
    ungroup()

  ### Heatmap dataframe for absence/presence of protein according to groups

  #  For each protein and group, determine if it's present (any non-NA)
  presence_matrix <- data_long %>%
    group_by(uniprot_id, group) %>%
    summarise(present = any(!is.na(Log2_LFQ_Intensity)), .groups = "drop") %>%
    dplyr::mutate(present = as.integer(present))  # 1 = present, 0 = all NA

  # From your presence_matrix
  proteins_absent_some_group <- presence_matrix %>%
    group_by(uniprot_id) %>%
    summarise(groups_absent = sum(present == 0), groups_present = sum(present == 1)) %>%
    filter(groups_absent >= 1, groups_present >= 1) %>%
    pull(uniprot_id)

  data_filtered <- data_long %>%
    filter(uniprot_id %in% proteins_absent_some_group)

  heatmap_matrix <- quant_data %>%
    as.data.frame() %>%
    filter(row.names(.) %in% unique(data_filtered[[col_pg]])) %>%
    rownames_to_column(col_pg) %>%
    left_join(prot_annot, by = col_pg) %>%
    dplyr::mutate(Unique_ID = paste("PG:", .data[[col_pg]], "Gene:", .data[[col_gene]])) %>%
    column_to_rownames("Unique_ID") %>%
    select(!colnames(prot_annot))

  list(
    quant_data = quant_data,
    meta_data = meta_data,
    prot_annot = prot_annot,
    data_long = data_long,
    vv_data = vv,
    stats_vv_data = stats_vv,
    rank_data = rank,
    corr_data = list(corr, p.mat),
    pca_data = p,
    heatmap_data = processed_matrix,
    vv_distrib_sample = valid_counts_all_complete,
    vv_distrib_grp = distribution_group_complete,
    heatmap_on_off = heatmap_matrix
  )


}

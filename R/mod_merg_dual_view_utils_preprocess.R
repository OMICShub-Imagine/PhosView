#' Preprocessing of Global + Phospho data for Dual View
#'
#' Combine global and phosphoproteomics statistics into a unified table
#' for the dual visualization dashboard.
#'
#' @param glob_stats List. The global analysis stats_data
#'   (glob$result()$stats_data)
#' @param phos_stats List. The phospho analysis stats_data
#'   (phos$result()$stats_data)
#'
#' @return A tibble/data.frame containing combined global + phospho data,
#'   with processed columns ready for dashboard display.
#'
#' @export
mod_merg_dual_view_utils_preprocess <- function(
    glob_stats,
    phos_stats,
    protein_id_col_glob,
    protein_id_col_phos,
    gene_col_glob = NULL,
    gene_col_phos = NULL,
    name_col_glob = NULL,
    name_col_phos = NULL,
    desc_col_glob = NULL,
    desc_col_phos = NULL,
    phospho_site_col
) {

  # ---------------- GLOBAL ----------------
  comparison_names_glob <- names(glob_stats$results)
  combined_glob <- setNames(
    lapply(comparison_names_glob, function(comp_name) glob_stats$results[[comp_name]]),
    comparison_names_glob
  )
  data_list_glob <- list(as.data.frame(glob_stats$data), glob_stats$annotation)
  global_data_final <- do.call(cbind, c(data_list_glob, combined_glob))
  colnames(global_data_final) <- paste(colnames(global_data_final), "GlobalAnalysis", sep = "_")
  rownames(global_data_final) <- NULL

  # ---------------- PHOSPHO ----------------
  comparison_names_phos <- names(phos_stats$results)
  combined_phos <- setNames(
    lapply(comparison_names_phos, function(comp_name) phos_stats$results[[comp_name]]),
    comparison_names_phos
  )
  data_list_phos <- list(phos_stats$data, phos_stats$annotation)
  phospho_data_final <- do.call(cbind, c(data_list_phos, combined_phos))
  colnames(phospho_data_final) <- paste(colnames(phospho_data_final), "PhosphoAnalysis", sep = "_")
  rownames(phospho_data_final) <- NULL


  # ---------------- MERGE ----------------
  data_final <- dplyr::full_join(
    x = global_data_final,
    y = phospho_data_final,
    by = setNames(
      paste0(protein_id_col_phos, "_PhosphoAnalysis"),
      paste0(protein_id_col_glob, "_GlobalAnalysis")
    ),
    keep = TRUE,
    multiple = "all"
  )

  # ---------------- COALESCE KEY PROTEIN COLUMNS ----------------
  if(!is.null(gene_col_glob) & !is.null(gene_col_phos)) {
    data_final <- data_final %>%
      dplyr::mutate(
        PG.Genes = dplyr::coalesce(
          !!rlang::sym(paste0(gene_col_glob, "_GlobalAnalysis")),
          !!rlang::sym(paste0(gene_col_phos, "_PhosphoAnalysis"))
        )
      )
  }

  if(!is.null(desc_col_glob) && length(desc_col_glob) > 0) {

    glob_desc_col <- paste0(desc_col_glob, "_GlobalAnalysis")

    phos_desc_col <- if(!is.null(desc_col_phos) && length(desc_col_phos) > 0) {
      paste0(desc_col_phos, "_PhosphoAnalysis")
    } else NULL

    if(!is.null(phos_desc_col) && phos_desc_col %in% colnames(data_final)) {
      # Phospho column exists → coalesce
      data_final <- data_final %>%
        dplyr::mutate(
          PG.ProteinDescriptions = dplyr::coalesce(
            !!rlang::sym(phos_desc_col),
            !!rlang::sym(glob_desc_col)
          )
        )
    } else {
      # Phospho column missing → just use global
      data_final <- data_final %>%
        dplyr::mutate(
          PG.ProteinDescriptions = !!rlang::sym(glob_desc_col)
        )
    }
  }

  if(!is.null(protein_id_col_glob) & !is.null(protein_id_col_phos)) {
    data_final <- data_final %>%
      dplyr::mutate(
        PG.ProteinGroups = dplyr::coalesce(
          !!rlang::sym(paste0(protein_id_col_glob, "_GlobalAnalysis")),
          !!rlang::sym(paste0(protein_id_col_phos, "_PhosphoAnalysis"))
        )
      )
  }

  if(!is.null(name_col_glob) & !is.null(name_col_phos)) {
    data_final <- data_final %>%
      dplyr::mutate(
        PG.ProteinNames = dplyr::coalesce(
          !!rlang::sym(paste0(name_col_glob, "_GlobalAnalysis")),
          !!rlang::sym(paste0(name_col_phos, "_PhosphoAnalysis"))
        )
      )
  }

  if(!is.null(phospho_site_col)) {
    data_final <- data_final %>%
      dplyr::mutate(
        PTM.CollapseKey = !!rlang::sym(paste0(phospho_site_col, "_PhosphoAnalysis"))
      )
  }

  # ---------------- SUBSET AND FORMAT ----------------
  subset_cols <- c("PTM.CollapseKey", "PG.Genes", "PG.ProteinGroups", "PG.ProteinDescriptions", "PG.ProteinNames")

  subset_data_final <- data_final %>%
    dplyr::select(
      dplyr::all_of(subset_cols),
      dplyr::contains("sig.FDR"),
      dplyr::contains("logFC"),
      dplyr::contains("P.Value"),
      dplyr::contains("adj.P.Val")
    )

  list(
    global_data_final  = global_data_final,
    phospho_data_final = phospho_data_final,
    subset_data_final  = subset_data_final
  )
}


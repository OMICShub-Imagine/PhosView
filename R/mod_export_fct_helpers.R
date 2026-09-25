## ============================================================
## Combine data + annotation
## ============================================================
combine_data_annotation <- function(x) {

  data_df <- as.data.frame(x$data) %>%
    tibble::rownames_to_column("uniprot_id")

  dplyr::left_join(x$annotation, data_df, by = "uniprot_id")
}

## ============================================================
## Combine stats_data (SUPPORTS MULTIPLE COMPARISONS)
## ============================================================
combine_stats <- function(stats_data) {

  # Annotation
  ann <- stats_data$annotation

  # Quantitative data
  quant <- as.data.frame(stats_data$data)

  # Statistics for each comparison
  stats_blocks <- lapply(names(stats_data$results), function(comp) {
    df <- stats_data$results[[comp]]
    # Prefix each column with the comparison name
    colnames(df) <- paste0(comp, "_", colnames(df))
    df
  })

  # Combine everything into one data.frame
  all_dfs <- c(list(ann, quant), stats_blocks)
  do.call(cbind, all_dfs)
}

## ============================================================
## Detect imputed cells
## ============================================================
get_imputed_matrix <- function(norm_mat, imput_mat) {
  is.na(norm_mat) & !is.na(imput_mat)
}

## ============================================================
## Clean annotation (REMOVE *_PTM_Group)
## ============================================================
clean_annotation <- function(annot) {
  dplyr::select(annot, -matches("_PTM_Group$"))
}

## ============================================================
## Combine data + cleaned annotation (phospho)
## ============================================================
combine_data_annotation_phos <- function(x) {

  data_df <- as.data.frame(x$data) |>
    tibble::rownames_to_column("uniprot_id")

  dplyr::left_join(
    clean_annotation(x$annotation),
    data_df,
    by = "uniprot_id"
  )
}

## ============================================================
## Combine stats_data (phospho, MULTIPLE COMPARISONS)
## ============================================================
combine_stats_phos <- function(stats_data) {

  # Raw quantitative data with rownames as column
  data_df <- as.data.frame(stats_data$data) |>
    tibble::rownames_to_column("uniprot_id")

  # Clean annotation
  annot_df <- clean_annotation(stats_data$annotation)

  # Merge annotation + data
  df <- dplyr::left_join(annot_df, data_df, by = "uniprot_id")

  # Add stats columns for each comparison
  if (!is.null(stats_data$results) && length(stats_data$results) > 0) {
    for (cmp in names(stats_data$results)) {
      cmp_df <- stats_data$results[[cmp]]

      # Prefix columns with comparison name
      colnames(cmp_df) <- paste(cmp, colnames(cmp_df), sep = "_")

      # Add to main df
      df <- cbind(df, cmp_df)
    }
  }

  df
}

## ============================================================
## Excel header styles
## ============================================================
.style_annotation <- openxlsx::createStyle(
  fgFill = "#BDD7EE",
  halign = "center",
  valign = "center",
  textDecoration = "bold"
)

.style_quant <- openxlsx::createStyle(
  fgFill = "#C6E0B4",
  halign = "center",
  valign = "center",
  textDecoration = "bold"
)

## ============================================================
## Stats styles (one color per comparison)
## ============================================================
make_stats_styles <- function(groups) {

  palette <- c(
    "#F8CBAD", "#D9D2E9", "#CFE2F3",
    "#FFF2CC", "#E2EFDA", "#FCE4D6"
  )

  if (length(groups) == 0) return(list())

  cols <- setNames(
    palette[seq_along(groups) %% length(palette) + 1],
    names(groups)
  )

  purrr::map(
    cols,
    ~ openxlsx::createStyle(
      fgFill = .x,
      halign = "center",
      valign = "center",
      textDecoration = "bold"
    )
  )
}

## ============================================================
## Add merged group headers
## ============================================================
add_grouped_headers <- function(
    wb, sheet, df,
    annot_cols,
    quant_cols,
    stats_groups = NULL
) {

  ## Annotation
  if (length(annot_cols)) {
    openxlsx::mergeCells(wb, sheet, cols = annot_cols, rows = 1)
    openxlsx::writeData(
      wb, sheet,
      "Protein annotation",
      startCol = min(annot_cols),
      startRow = 1
    )
    openxlsx::addStyle(wb, sheet, .style_annotation, 1, annot_cols, TRUE)
  }

  ## Quantitative
  if (length(quant_cols)) {
    openxlsx::mergeCells(wb, sheet, cols = quant_cols, rows = 1)
    openxlsx::writeData(
      wb, sheet,
      "Quantitative data for each sample",
      startCol = min(quant_cols),
      startRow = 1
    )
    openxlsx::addStyle(wb, sheet, .style_quant, 1, quant_cols, TRUE)
  }

  ## Statistics
  if (!is.null(stats_groups)) {

    stats_styles <- make_stats_styles(stats_groups)

    for (grp in names(stats_groups)) {

      cols <- stats_groups[[grp]]

      openxlsx::mergeCells(wb, sheet, cols = cols, rows = 1)
      openxlsx::writeData(
        wb, sheet,
        grp,
        startCol = min(cols),
        startRow = 1
      )
      openxlsx::addStyle(
        wb, sheet,
        stats_styles[[grp]],
        rows = 1,
        cols = cols,
        gridExpand = TRUE
      )
    }
  }
}
## ============================================================
## Column detection helpers
## ============================================================
get_annot_cols <- function(x, phospho = FALSE) {
  annot <- if (phospho) clean_annotation(x$annotation) else x$annotation
  seq_len(ncol(annot))
}

get_quant_cols <- function(df, x) {
  cols <- which(colnames(df) %in% colnames(x$data))

  if (length(cols) > 1) {
    cols <- min(cols):max(cols)  # force contiguous range
  }

  cols
}

get_stats_groups <- function(df) {

  stat_names <- c(
    "logFC", "CI.L", "CI.R", "average_intensity",
    "t", "B", "P.Value", "adj.P.Val", "sig.PVal", "sig.FDR"
  )

  stat_pattern <- paste0("_(", paste(stat_names, collapse = "|"), ")$")

  stat_cols <- which(grepl(stat_pattern, colnames(df)))

  split(
    stat_cols,
    sub(stat_pattern, "", colnames(df)[stat_cols])
  )
}

extract_ora_tables <- function(ora_obj) {
  # ora_obj: result of run_all_ORA() for one mode
  out <- list()

  for (cmp in names(ora_obj)) {
    for (db in names(ora_obj[[cmp]])) {

      enr <- ora_obj[[cmp]][[db]]

      if (inherits(enr, "enrichResult")) {
        out[[paste(cmp, db, sep = "_")]] <- enr@result
      }
    }
  }

  out
}

extract_signalome_tables <- function(signalome) {

  out <- list()

  ## -------------------------------------------------
  ## Utility: generate Excel-safe sheet names (≤31 chars)
  ## -------------------------------------------------
  make_sheet_name <- function(cmp, direction = NULL, type) {

    # Shorten comparison names
    cmp <- gsub("_vs_", "v", cmp)   # PL_vs_CTRL → PLvCTRL

    # Shorten content types
    type <- switch(
      type,
      "UP_filt_phospho"     = "UP_raw",
      "DOWN_filt_phospho"   = "DOWN_raw",
      "UP_standardized"     = "UP_z",
      "DOWN_standardized"   = "DOWN_z",
      "ksActivityMatrix"    = "KS_act",
      "combinedScoreMatrix" = "KS_comb",
      "UP_kin_pred"         = "UP_pred",
      "DOWN_kin_pred"       = "DOWN_pred",
      type
    )

    nm <- paste(c(cmp, direction, type), collapse = "_")

    substr(nm, 1, 31)
  }

  ## -------------------------------------------------
  ## Utility: extract score matrices safely
  ## -------------------------------------------------
  extract_score_tables <- function(lst, direction, score_name) {

    res <- list()

    for (cmp in names(lst)) {

      obj <- lst[[cmp]]

      if (!is.null(obj[[score_name]]) && nrow(obj[[score_name]]) > 0) {

        sheet <- make_sheet_name(cmp, direction, score_name)
        res[[sheet]] <- as.data.frame(obj[[score_name]])
      }
    }

    res
  }

  ## -------------------------------------------------
  ## 1. Filtered phosphosites (raw intensities)
  ## -------------------------------------------------
  for (cmp in names(signalome$positive_filtered_data)) {
    out[[make_sheet_name(cmp, NULL, "UP_filt_phospho")]] <-
      as.data.frame(signalome$positive_filtered_data[[cmp]])
  }

  for (cmp in names(signalome$negative_filtered_data)) {
    out[[make_sheet_name(cmp, NULL, "DOWN_filt_phospho")]] <-
      as.data.frame(signalome$negative_filtered_data[[cmp]])
  }

  ## -------------------------------------------------
  ## 2. Standardized phosphosite matrices
  ## -------------------------------------------------
  for (cmp in names(signalome$positive_standardized_data)) {
    out[[make_sheet_name(cmp, NULL, "UP_standardized")]] <-
      as.data.frame(signalome$positive_standardized_data[[cmp]])
  }

  for (cmp in names(signalome$negative_standardized_data)) {
    out[[make_sheet_name(cmp, NULL, "DOWN_standardized")]] <-
      as.data.frame(signalome$negative_standardized_data[[cmp]])
  }

  ## -------------------------------------------------
  ## 3. Kinase scores (FIXED)
  ## -------------------------------------------------
  out <- c(
    out,
    extract_score_tables(
      signalome$kinase_score_mat_list_UP,
      "UP",
      "ksActivityMatrix"
    ),
    extract_score_tables(
      signalome$kinase_score_mat_list_DOWN,
      "DOWN",
      "ksActivityMatrix"
    ),
    extract_score_tables(
      signalome$kinase_score_mat_list_UP,
      "UP",
      "combinedScoreMatrix"
    ),
    extract_score_tables(
      signalome$kinase_score_mat_list_DOWN,
      "DOWN",
      "combinedScoreMatrix"
    )
  )

  ## -------------------------------------------------
  ## 4. Kinase–substrate prediction matrices
  ## -------------------------------------------------
  for (cmp in names(signalome$kinase_predmat_list_UP)) {

    pred <- signalome$kinase_predmat_list_UP[[cmp]]

    if (!is.null(pred) && nrow(pred) > 0) {
      out[[make_sheet_name(cmp, "UP", "UP_kin_pred")]] <- as.data.frame(pred)
    }
  }

  for (cmp in names(signalome$kinase_predmat_list_DOWN)) {

    pred <- signalome$kinase_predmat_list_DOWN[[cmp]]

    if (!is.null(pred) && nrow(pred) > 0) {
      out[[make_sheet_name(cmp, "DOWN", "DOWN_kin_pred")]] <- as.data.frame(pred)
    }
  }

  out
}




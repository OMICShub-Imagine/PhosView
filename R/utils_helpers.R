harmonize_annotation <- function(df, type = c("glob", "phos")) {

  type <- match.arg(type)

  if (type == "glob") {

    if ("Genes" %in% names(df)) {
      df$gene <- df$Genes
    } else if ("PG.Genes" %in% names(df)) {
      df$gene <- df$`PG.Genes`
    }

    if ("Protein.Group" %in% names(df)) {
      df$protein_group <- df$`Protein.Group`
    } else if ("PG.ProteinGroups" %in% names(df)) {
      df$protein_group <- df$`PG.ProteinGroups`
    }

    if ("Protein.Names" %in% names(df)) {
      df$protein_name <- df$`Protein.Names`
    } else if ("PG.ProteinNames" %in% names(df)) {
      df$protein_name <- df$`PG.ProteinNames`
    }

    if ("First.Protein.Description" %in% names(df)) {
      df$description <- df$`First.Protein.Description`
    } else if ("PG.ProteinDescriptions" %in% names(df)) {
      df$description <- df$`PG.ProteinDescriptions`
    }

  }

  if (type == "phos") {

    if ("Gene.Names" %in% names(df)) {
      df$gene <- df$`Gene.Names`
    } else if ("PG.Genes" %in% names(df)) {
      df$gene <- df$`PG.Genes`
    }

    #  unify site_id → phosphosite
    if ("PTM.CollapseKey" %in% names(df)) {
      df$phosphosite <- df$`PTM.CollapseKey`
    } else if ("uniprot_id" %in% names(df)) {
      df$phosphosite <- df$uniprot_id
    }

    #  unify site_id → phosphosite
    if ("PTM.FlankingRegion" %in% names(df)) {
      df$FlankingRegionSeq <- df$`PTM.FlankingRegion`
    } else if ("Sequence" %in% names(df)) {
      df$FlankingRegionSeq <- df$Sequence
    }

  }

  df
}


safe_enrich_df <- function(x) {
  if (is.null(x) || is.null(x@result) || nrow(x@result) == 0) {
    return(data.frame(Message = "No enrichment found"))
  }
  as.data.frame(x@result)
}

has_sig <- function(res) {
  if (is.null(res) || is.null(res@result)) return(FALSE)
  any(res@result$p.adjust < 0.05, na.rm = TRUE)
}

safe_plot <- function(res, plot_fun) {
  if (!has_sig(res)) {
    plot.new()
    text(0.5, 0.5, "No significant results to display")
  } else {
    plot_fun(res)
  }
}

make_termsim <- function(res) {
  if (is.null(res) || is.null(res@result) || nrow(res@result) < 2) return(NULL)
  tryCatch(enrichplot::pairwise_termsim(res), error = function(e) NULL)
}

# Adapted from proteoDA
# Original source:
# https://github.com/ByrumLab/proteoDA
# Licensed under the GNU General Public License v3.0
validate_DAList_local <- function(x) {
  `%||%` <- function(a, b) if (is.null(a)) b else a

  # ---- Required / optional slots ----
  required_slots <- c("data", "annotation", "metadata", "design", "eBayes_fit", "results", "tags")
  optional_slots <- c("filtered_proteins_per_contrast", "data_per_contrast", "annotation_per_contrast")

  missing_required <- setdiff(required_slots, names(x))
  if (length(missing_required) > 0) {
    cli::cli_abort("The DAList is missing the following required slot{?s}: {missing_required}")
  }

  # ---- Basic shape checks ----
  if (any(c("tbl_df", "tbl") %in% class(x$data)))
    cli::cli_abort("Data slot is a tibble, must be a matrix or data frame")
  if (any(c("tbl_df", "tbl") %in% class(x$annotation)))
    cli::cli_abort("Annotation slot is a tibble, must be a matrix or data frame")
  if (any(c("tbl_df", "tbl") %in% class(x$metadata)))
    cli::cli_abort("Metadata slot is a tibble, must be a matrix or data frame")

  if (!any(c(is.data.frame(x$data), is.matrix(x$data))))
    cli::cli_abort("The 'data' slot must be a data frame or matrix")
  if (!all(apply(x$data, 2, is.numeric)))
    cli::cli_abort("The 'data' slot must contain only numeric data")

  if (!"uniprot_id" %in% colnames(x$annotation))
    cli::cli_abort("The annotation data must contain a column named 'uniprot_id'")
  if (any(duplicated(x$annotation$uniprot_id)))
    cli::cli_abort("The 'uniprot_id' values in annotation must be unique")
  if (nrow(x$data) != nrow(x$annotation))
    cli::cli_abort("'data' and 'annotation' must have the same number of rows")
  if (!identical(rownames(x$data), rownames(x$annotation)))
    cli::cli_abort("Rownames for 'data' and 'annotation' must match")
  if (!identical(rownames(x$data), x$annotation$uniprot_id))
    cli::cli_abort("Rownames for 'data'/'annotation' must equal the 'uniprot_id' column")

  # metadata to data samples
  if (nrow(x$metadata) != ncol(x$data)) {
    cli::cli_abort("Sample count mismatch: metadata rows = {nrow(x$metadata)}, data columns = {ncol(x$data)}")
  }
  if (!identical(colnames(x$data), rownames(x$metadata))) {
    miss_in_md  <- setdiff(colnames(x$data), rownames(x$metadata))
    miss_in_dat <- setdiff(rownames(x$metadata), colnames(x$data))
    same_set    <- setequal(colnames(x$data), rownames(x$metadata))
    hint <- if (same_set) " (same samples, different order -- reorder metadata to match data columns)" else ""
    det <- c(
      if (length(miss_in_md))  paste0("Missing in metadata: ", paste(miss_in_md, collapse = ", ")),
      if (length(miss_in_dat)) paste0("Missing in data: ",     paste(miss_in_dat, collapse = ", "))
    )
    cli::cli_abort(c(
      "Row/column name mismatch between data and metadata{hint}",
      if (length(det)) setNames(det, rep("x", length(det)))
    ))
  }

  # ---- Design checks ----
  if (!is.null(x$design)) {
    need <- c("design_matrix", "design_formula")
    if (length(x$design) < 2)
      cli::cli_abort("'design' must be a list of at least 2 items: design_formula and design_matrix")
    if (!all(need %in% names(x$design)))
      cli::cli_abort("'design' must include 'design_matrix' and 'design_formula'")

    # assign: warn (allow custom matrices)
    if (!"assign" %in% names(attributes(x$design$design_matrix))) {
      cli::cli_warn("'design_matrix' has no 'assign' attribute; recommended to build via model.matrix.")
    }

    # design_matrix to metadata alignment
    X  <- x$design$design_matrix
    md <- x$metadata
    if (nrow(X) != nrow(md)) {
      cli::cli_abort("Row count mismatch: design_matrix rows = {nrow(X)}, metadata rows = {nrow(md)}")
    }
    if (!identical(rownames(X), rownames(md))) {
      miss_in_X  <- setdiff(rownames(md), rownames(X))
      miss_in_md <- setdiff(rownames(X), rownames(md))
      same_set   <- setequal(rownames(X), rownames(md))
      hint <- if (same_set) " (same samples, different order -- reorder design_matrix to match metadata)" else ""
      det <- c(
        if (length(miss_in_X))  paste0("Missing in design_matrix: ", paste(miss_in_X,  collapse = ", ")),
        if (length(miss_in_md)) paste0("Missing in metadata: ",      paste(miss_in_md, collapse = ", "))
      )
      cli::cli_abort(c(
        "Row names of design_matrix must match row names of metadata{hint}",
        if (length(det)) setNames(det, rep("x", length(det)))
      ))
    }

    # # random factor present in metadata
    # if (!is.null(x$design$random_factor)) {
    #   if (!x$design$random_factor %in% colnames(x$metadata))
    #     cli::cli_abort("Random factor '{x$design$random_factor}' not found in metadata")
    # }
    # random factor present in metadata
    if (!is.null(x$design$random_factor)) {
      if (!x$design$random_factor %in% colnames(x$metadata)) {
        cli::cli_abort(
          "Design includes random factor term '{x$design$random_factor}', but metadata has no column with that name"
        )
      }
    }


    # contrast info: require both or neither
    if (sum(c("contrast_matrix", "contrast_vector") %in% names(x$design)) == 1) {
      cli::cli_abort("If contrast information is included, you must have  both 'contrast_matrix' and 'contrast_vector' slots")
    }

    # contrast_matrix to design_matrix columns
    if ("contrast_matrix" %in% names(x$design)) {
      cm <- x$design$contrast_matrix
      if (!identical(rownames(cm), colnames(X))) {
        miss_in_cm <- setdiff(colnames(X), rownames(cm))
        miss_in_X  <- setdiff(rownames(cm), colnames(X))
        same_set   <- setequal(rownames(cm), colnames(X))
        hint <- if (same_set) " (same terms, different order -- reorder contrast_matrix rows to match design_matrix columns)" else ""
        det <- c(
          if (length(miss_in_cm)) paste0("Missing in contrast_matrix rows: ", paste(miss_in_cm, collapse = ", ")),
          if (length(miss_in_X))  paste0("Missing in design_matrix cols: ",  paste(miss_in_X,  collapse = ", "))
        )
        cli::cli_abort(c(
          "contrast_matrix rows do not match design_matrix columns{hint}",
          if (length(det)) setNames(det, rep("x", length(det)))
        ))
      }
    }
  }

  # ---- eBayes fit checks ----
  if (!is.null(x$eBayes_fit)) {

    # Coerce to a list of fits and do type checks
    if (inherits(x$eBayes_fit, "MArrayLM")) {
      fits <- list(x$eBayes_fit)
    } else if (is.list(x$eBayes_fit)) {
      ok <- vapply(
        x$eBayes_fit,
        function(obj) inherits(obj, "MArrayLM"),
        FUN.VALUE = logical(1)
      )
      if (!all(ok)) {
        cli::cli_abort("eBayes_fit must be of class 'MArrayLM' or a list of 'MArrayLM' objects")
      }
      fits <- x$eBayes_fit
    } else {
      cli::cli_abort("eBayes_fit must be of class 'MArrayLM' or a list of 'MArrayLM' objects")
    }

    # ---- Random factor consistency checks ----
    # design-side declaration
    rf_design <- x$design$random_factor %||% NULL

    # fit-side declaration (set in fit_limma_model via attr(efit, "random_factor_name"))
    rf_fit_name <- attr(fits[[1]], "random_factor_name", exact = TRUE)

    # Case 2: design claims a random factor, but model fit doesn't have one
    if (!is.null(rf_design) && is.null(rf_fit_name)) {
      cli::cli_abort(
        "Design includes a random factor term '{rf_design}', but the model fit does not include a corresponding random-effects structure.",
        rf_design = rf_design
      )
    }

    # Case 3: model fit has a random factor, but design tag is missing
    if (is.null(rf_design) && !is.null(rf_fit_name)) {
      cli::cli_abort(
        "Model fit includes a random-effects structure (random factor '{rf_fit_name}'), but 'design$random_factor' is not set.",
        rf_fit_name = rf_fit_name
      )
    }
  }


  # ---- Results checks (lenient on labels) ----
  if (!is.null(x$results)) {
    term_names <- if (!is.null(x$design$contrast_matrix)) {
      colnames(x$design$contrast_matrix)
    } else {
      colnames(x$design$design_matrix)
    }

    expected <- if (!is.null(x$tags$extract_intercept) && isTRUE(x$tags$extract_intercept)) {
      term_names
    } else {
      term_names[!grepl("ntercept", term_names, ignore.case = TRUE)]
    }

    res_names <- names(x$results)

    if (!identical(res_names, expected)) {
      unmapped <- setdiff(res_names, expected)
      if (length(unmapped) > 0) {
        dcols <- colnames(x$design$design_matrix)
        still_bad <- character(0)
        for (lab in unmapped) {
          ci <- tryCatch(x$tags$per_contrast[[lab]]$contrast_info, error = function(e) NULL)
          ok <- FALSE
          if (!is.null(ci) && !is.null(ci$contrast_expression)) {
            expr <- as.character(ci$contrast_expression)[1]
            tokens <- unique(unlist(regmatches(expr, gregexpr("[A-Za-z0-9_.:]+", expr))))
            if (length(intersect(tokens, dcols)) > 0) ok <- TRUE
          }
          if (!ok) still_bad <- c(still_bad, lab)
        }
        if (length(still_bad) > 0) {
          cli::cli_warn("Results names do not match design/contrast terms and lack usable tags: {still_bad}. Proceeding.")
        }
      } else {
        #cli::cli_warn("Results names differ from expected design/contrast terms; proceeding with arbitrary labels.")
        cli::cli_inform("Results names differ from expected design/contrast terms; proceeding with arbitrary labels.")
      }
    }

    # Per-contrast row checks: must be subset of data rows (order may differ)
    for (term in res_names) {
      result_rows <- rownames(x$results[[term]])
      data_rows   <- rownames(x$data)

      if (!all(result_rows %in% data_rows)) {
        cli::cli_abort("Some rows in results for contrast '{term}' are not present in 'data'", term = term)
      }

      if ("filtered_proteins_per_contrast" %in% names(x) && !is.null(x$filtered_proteins_per_contrast[[term]])) {
        expected_subset <- x$filtered_proteins_per_contrast[[term]]
        if (!all(result_rows %in% expected_subset)) {
          cli::cli_warn("Some rows in results for contrast '{term}' do not match the filtered protein list", term = term)
        }
      }
    }
  }

  # ---- Normalize order & class (back-compat for callers expecting an object) ----
  ordered_names <- c(required_slots, intersect(optional_slots, names(x)))
  ordered <- x[ordered_names]
  extras  <- setdiff(names(x), ordered_names)
  if (length(extras)) ordered <- c(ordered, x[extras])  # preserve extra/custom slots

  class(ordered) <- unique(c("DAList", class(ordered)))
  ordered
}

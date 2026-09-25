#' Compute signalome and kinase–substrate inference from phosphoproteomics data
#'
#' Processes phosphoproteomics differential analysis results to identify
#' up- and down-regulated phosphosites, standardize signals, extract PTM
#' sequences, infer kinase–substrate scores, and predict kinase activities
#' using \pkg{PhosR}.
#'
#' @param phospho_stats A list containing phosphoproteomics differential
#'   analysis results. Must include \code{$results}, \code{$data}, and
#'   \code{$annotation}.
#' @param species Character string specifying the organism
#'   (e.g. \code{"human"}, \code{"mouse"}).
#' @param PhosphoSite.Db A kinase–substrate database compatible with
#'   \code{\link[PhosR]{kinaseSubstrateScore}}.
#'
#' @details
#' The function performs the following steps for each comparison:
#' \enumerate{
#'   \item Filters significantly up- and down-regulated phosphosites
#'   \item Standardizes phosphosite intensities
#'   \item Extracts PTM flanking sequences
#'   \item Computes kinase–substrate score matrices
#'   \item Predicts kinase–substrate interactions
#' }
#'
#' All failures in kinase scoring or prediction are safely handled and return
#' empty result structures instead of stopping execution.
#'
#' @return A list containing:
#' \describe{
#'   \item{positive_filtered_data}{Filtered up-regulated phosphosites}
#'   \item{negative_filtered_data}{Filtered down-regulated phosphosites}
#'   \item{positive_standardized_data}{Standardized up-regulated matrices}
#'   \item{negative_standardized_data}{Standardized down-regulated matrices}
#'   \item{positive_PTM_seq_list}{PTM sequences for up-regulated sites}
#'   \item{negative_PTM_seq_list}{PTM sequences for down-regulated sites}
#'   \item{kinase_score_mat_list_UP}{Kinase–substrate score results (UP)}
#'   \item{kinase_score_mat_list_DOWN}{Kinase–substrate score results (DOWN)}
#'   \item{kinase_predmat_list_UP}{Predicted kinase substrates (UP)}
#'   \item{kinase_predmat_list_DOWN}{Predicted kinase substrates (DOWN)}
#' }
#'
#' @seealso
#' \code{\link[PhosR]{standardise}},
#' \code{\link[PhosR]{kinaseSubstrateScore}},
#' \code{\link[PhosR]{kinaseSubstratePred}}
#'
#'@importFrom PhosR standardise kinaseSubstrateScore kinaseSubstratePred
#'
#' @export
compute_signalome <- function(phospho_stats, species, PhosphoSite.Db) {

  ### ------------------------------------------------------------
  ### 1. FILTER PHOSPHOSITES (UP / DOWN)
  ### ------------------------------------------------------------

  positive_filtered_data <- list()
  negative_filtered_data <- list()

  for (comparison_name in names(phospho_stats$results)) {

    comp <- phospho_stats$results[[comparison_name]]

    up_phosphosite   <- rownames(comp %>% filter(sig.FDR == 1))
    down_phosphosite <- rownames(comp %>% filter(sig.FDR == -1))

    positive_filtered_data[[comparison_name]] <-
      phospho_stats$data[rownames(phospho_stats$data) %in% up_phosphosite, ]

    negative_filtered_data[[comparison_name]] <-
      phospho_stats$data[rownames(phospho_stats$data) %in% down_phosphosite, ]
  }

  ### ------------------------------------------------------------
  ### 2. STANDARDISE
  ### ------------------------------------------------------------

  positive_standardized_data <- list()
  negative_standardized_data <- list()

  for (comparison_name in names(positive_filtered_data)) {
    df <- positive_filtered_data[[comparison_name]]
    positive_standardized_data[[comparison_name]] <- PhosR::standardise(as.matrix(df))
  }

  for (comparison_name in names(negative_filtered_data)) {
    df <- negative_filtered_data[[comparison_name]]
    negative_standardized_data[[comparison_name]] <- PhosR::standardise(as.matrix(df))
  }

  ### ------------------------------------------------------------
  ### 3. EXTRACT PTM SEQUENCES
  ### ------------------------------------------------------------

  extract_PTM_seq <- function(rows) {
    phospho_stats$annotation[rownames(phospho_stats$annotation) %in% rows,
                             "FlankingRegionSeq", drop = FALSE]
  }

  positive_PTM_seq_list <- list()
  negative_PTM_seq_list <- list()

  for (comparison_name in names(positive_standardized_data)) {
    rn <- rownames(positive_standardized_data[[comparison_name]])
    positive_PTM_seq_list[[comparison_name]] <- extract_PTM_seq(rn)
  }

  for (comparison_name in names(negative_standardized_data)) {
    rn <- rownames(negative_standardized_data[[comparison_name]])
    negative_PTM_seq_list[[comparison_name]] <- extract_PTM_seq(rn)
  }

  ### ------------------------------------------------------------
  ### 4. MODIFY ROWNAMES
  ### ------------------------------------------------------------

  modify_row_names <- function(row_names) {
    sapply(strsplit(row_names, "_"), function(x) {
      if (grepl("^M\\d+$", tail(x, 1))) {
        x <- head(x, -1)
      }
      paste0(paste(x, collapse = ";"), ";")
    })
  }

  for (comparison_name in names(positive_standardized_data)) {
    mat <- positive_standardized_data[[comparison_name]]
    if (nrow(mat) > 0) {
      rownames(mat) <- modify_row_names(rownames(mat))
    }
    positive_standardized_data[[comparison_name]] <- mat
  }

  for (comparison_name in names(negative_standardized_data)) {
    mat <- negative_standardized_data[[comparison_name]]
    if (nrow(mat) > 0) {
      rownames(mat) <- modify_row_names(rownames(mat))
    }
    negative_standardized_data[[comparison_name]] <- mat
  }

  ### ------------------------------------------------------------
  ### 5. KINASE SUBSTRATE SCORE (UP / DOWN)
  ### ------------------------------------------------------------

  kinase_score_mat_list_UP <- list()
  kinase_score_mat_list_DOWN <- list()

  for (comparison_name in names(positive_standardized_data)) {

    mat <- positive_standardized_data[[comparison_name]]
    seqs <- as.list(positive_PTM_seq_list[[comparison_name]]$FlankingRegionSeq)

    # uppercase gene names
    rownames(mat) <- toupper(rownames(mat))

    res <- tryCatch({
      PhosR::kinaseSubstrateScore(
        substrate.list = PhosphoSite.Db,
        mat = mat,
        seqs = seqs,
        numMotif = 5,
        numSub = 1,
        species = species,
        verbose = FALSE
      )
    }, error = function(e) {
      list(
        motifScoreMatrix = data.frame(),
        profileScoreMatrix = data.frame(),
        combinedScoreMatrix = data.frame(),
        ksActivityMatrix = data.frame(),
        weights = data.frame()
      )
    })

    kinase_score_mat_list_UP[[comparison_name]] <- res
  }

  for (comparison_name in names(negative_standardized_data)) {

    mat <- negative_standardized_data[[comparison_name]]
    seqs <- as.list(negative_PTM_seq_list[[comparison_name]]$FlankingRegionSeq)

    rownames(mat) <- toupper(rownames(mat))

    res <- tryCatch({

      # Run in a temp dir, automatically restores WD after execution
      res <- withr::with_tempdir({
        suppressWarnings(
          PhosR::kinaseSubstrateScore(
            substrate.list = PhosphoSite.Db,
            mat = mat,
            seqs = seqs,
            numMotif = 5,
            numSub = 1,
            species = species,
            verbose = FALSE
          )
        )
      })

      res

    }, error = function(e) {
      list(
        motifScoreMatrix = data.frame(),
        profileScoreMatrix = data.frame(),
        combinedScoreMatrix = data.frame(),
        ksActivityMatrix = data.frame(),
        weights = data.frame()
      )
    })


    kinase_score_mat_list_DOWN[[comparison_name]] <- res
  }

  ### ------------------------------------------------------------
  ### 6. KINASE PREDICTIONS (UP / DOWN)
  ### ------------------------------------------------------------

  set.seed(1)
  kinase_predmat_list_UP <- list()
  kinase_predmat_list_DOWN <- list()

  for (comparison_name in names(kinase_score_mat_list_UP)) {
    kmat <- kinase_score_mat_list_UP[[comparison_name]]

    pred <- tryCatch({
      PhosR::kinaseSubstratePred(
        kmat,
        ensembleSize = 10,
        top = 30,
        cs = 0.8,
        inclusion = 5,
        iter = 5,
        verbose = FALSE
      )
    }, error = function(e) NULL)

    kinase_predmat_list_UP[[comparison_name]] <- pred
  }

  for (comparison_name in names(kinase_score_mat_list_DOWN)) {
    kmat <- kinase_score_mat_list_DOWN[[comparison_name]]

    pred <- tryCatch({
      PhosR::kinaseSubstratePred(
        kmat,
        ensembleSize = 10,
        top = 30,
        cs = 0.8,
        inclusion = 5,
        iter = 5,
        verbose = FALSE
      )
    }, error = function(e) NULL)

    kinase_predmat_list_DOWN[[comparison_name]] <- pred
  }

  ### ------------------------------------------------------------
  ### RETURN ALL OBJECTS AS A SINGLE LIST
  ### ------------------------------------------------------------

  return(list(
    positive_filtered_data = positive_filtered_data,
    negative_filtered_data = negative_filtered_data,
    positive_standardized_data = positive_standardized_data,
    negative_standardized_data = negative_standardized_data,
    positive_PTM_seq_list = positive_PTM_seq_list,
    negative_PTM_seq_list = negative_PTM_seq_list,
    kinase_score_mat_list_UP = kinase_score_mat_list_UP,
    kinase_score_mat_list_DOWN = kinase_score_mat_list_DOWN,
    kinase_predmat_list_UP = kinase_predmat_list_UP,
    kinase_predmat_list_DOWN = kinase_predmat_list_DOWN
  ))
}


#' Render interactive kinase–substrate heatmap using heatmaply
#'
#' Creates an interactive heatmap of kinase–substrate combined scores,
#' displaying the top phosphosites per kinase and annotating kinases
#' by group and family.
#'
#' @param phosScoringMatrices A list containing kinase–substrate scoring
#'   results, including \code{combinedScoreMatrix}.
#' @param KinaseFamily A data frame describing kinase families and groups.
#'   Row names must match kinase names.
#' @param comparison Character string identifying the comparison being plotted.
#'
#' @details
#' The function selects the top three phosphosites per kinase based on
#' combined scores. If no valid data are available, a placeholder heatmap
#' is returned.
#'
#' @return An interactive \code{heatmaply} object.
#'
#' @seealso \code{\link[heatmaply]{heatmaply}}
#'
#' @export

renderHeatmap_heatmaply <- function(phosScoringMatrices, KinaseFamily, comparison) {

  #------------------------------------------------------------
  # Helper
  #------------------------------------------------------------
  empty_heatmap <- function(msg = "No data") {
    heatmaply::heatmaply(
      matrix(0, 1, 1),
      row_names = FALSE,
      col_names = FALSE,
      main = msg
    )
  }

  #------------------------------------------------------------
  # Check input
  #------------------------------------------------------------
  if (is.null(phosScoringMatrices) ||
      is.null(phosScoringMatrices$combinedScoreMatrix)) {

    message("combinedScoreMatrix not found for comparison ", comparison)
    return(empty_heatmap("No data"))
  }

  mat <- phosScoringMatrices$combinedScoreMatrix

  if (!is.matrix(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    message("Empty score matrix for ", comparison)
    return(empty_heatmap("No data"))
  }

  #------------------------------------------------------------
  # Extract top 3 sites per kinase
  #------------------------------------------------------------
  sites <- character(0)

  for (i in seq_len(ncol(mat))) {

    top_n <- min(3, nrow(mat))

    top3 <- names(
      sort(mat[, i], decreasing = TRUE)[seq_len(top_n)]
    )

    sites <- union(sites, top3)
  }

  #------------------------------------------------------------
  # Keep kinases present in KinaseFamily
  #------------------------------------------------------------
  o <- intersect(colnames(mat), rownames(KinaseFamily))

  if (length(sites) == 0 || length(o) == 0) {
    message("No valid sites or kinases for heatmap in ", comparison)
    return(empty_heatmap("No data"))
  }

  submat <- mat[sites, o, drop = FALSE]

  message(
    "Heatmap dimensions for ",
    comparison,
    ": ",
    nrow(submat),
    " rows x ",
    ncol(submat),
    " columns"
  )

  #------------------------------------------------------------
  # Not enough objects for clustering
  #------------------------------------------------------------
  row_cluster <- nrow(submat) >= 2
  col_cluster <- ncol(submat) >= 2

  if (!row_cluster && !col_cluster) {

    message(
      "Not enough rows/columns for clustering in ",
      comparison
    )

    return(
      heatmaply::heatmaply(
        submat,
        Rowv = FALSE,
        Colv = FALSE,
        main = paste("Insufficient data for clustering -", comparison)
      )
    )
  }

  #------------------------------------------------------------
  # Annotation
  #------------------------------------------------------------
  group_vec <- as.character(KinaseFamily[o, "kinase_group"])
  family_vec <- as.character(KinaseFamily[o, "kinase_family"])

  group_vec[is.na(group_vec)] <- "Unknown"
  family_vec[is.na(family_vec)] <- "Unknown"

  side_ann <- data.frame(
    Kinase_Group  = factor(group_vec,
                           levels = names(KINASE_GROUP_COLORS)),
    Kinase_Family = factor(family_vec,
                           levels = names(KINASE_FAMILY_COLORS)),
    row.names = o,
    stringsAsFactors = FALSE
  )

  #------------------------------------------------------------
  # Render heatmap
  #------------------------------------------------------------
  heatmaply::heatmaply(
    submat,
    scale = "none",
    Rowv = row_cluster,
    Colv = col_cluster,
    colors = grDevices::colorRampPalette(
      c("blue", "yellow", "red")
    )(200),
    main = "Top 3 phosphosite(s) for each kinase",
    row_side_colors = NULL,
    col_side_colors = side_ann,
    show_dendrogram = c(row_cluster, col_cluster),
    hide_colorbar = FALSE
  )
}

#' Render kinase–substrate heatmap for up-regulated phosphosites
#'
#' Convenience wrapper around \code{renderHeatmap_heatmaply()} to visualize
#' kinase–substrate scores derived from up-regulated phosphosites.
#'
#' @param kinsub_results Output list from \code{compute_signalome()}.
#' @param comparison Character string identifying the comparison to plot.
#'
#' @return An interactive \code{heatmaply} object.
#'
#' @export
renderHeatmap_UP <- function(kinsub_results, comparison) {
  phosScoringMatrices <- kinsub_results$kinase_score_mat_list_UP[[comparison]]
  renderHeatmap_heatmaply(phosScoringMatrices, KinaseFamily, comparison)
}

#' Render kinase–substrate heatmap for down-regulated phosphosites
#'
#' Convenience wrapper around \code{renderHeatmap_heatmaply()} to visualize
#' kinase–substrate scores derived from down-regulated phosphosites.
#'
#' @param kinsub_results Output list from \code{compute_signalome()}.
#' @param comparison Character string identifying the comparison to plot.
#'
#' @return An interactive \code{heatmaply} object.
#'
#' @export
renderHeatmap_DOWN <- function(kinsub_results, comparison) {
  phosScoringMatrices <- kinsub_results$kinase_score_mat_list_DOWN[[comparison]]
  renderHeatmap_heatmaply(phosScoringMatrices, KinaseFamily, comparison)
}

# renderHeatmap_UP <- function(kinsub_results, comparison) {
#
#   # Read the RDS file reactively
#   phosScoringMatrices <- kinsub_results$kinase_score_mat_list_UP[[comparison]]
#
#   # If phosScoringMatrices is NULL, return a placeholder Heatmap
#   if (is.null(phosScoringMatrices) || is.null(phosScoringMatrices$combinedScoreMatrix)) {
#     message("combinedScoreMatrix not found in the data for comparison ", comparison)
#     return(ComplexHeatmap::Heatmap(matrix(0, nrow = 1, ncol = 1), name = "No Data"))
#   }
#
#   # Extract top sites
#   sites <- c()
#   for (i in seq_len(ncol(phosScoringMatrices$combinedScoreMatrix))) {
#     sites <- union(sites, names(
#       sort(phosScoringMatrices$combinedScoreMatrix[, i], decreasing = TRUE)[seq_len(3)]
#     ))
#   }
#
#   o <- intersect(colnames(phosScoringMatrices$combinedScoreMatrix),
#                  rownames(KinaseFamily))
#
#   # Get unique values for group and family
#   unique_groups <- unique(KinaseFamily[o, "kinase_group"])
#   unique_families <- unique(KinaseFamily[o, "kinase_family"])
#
#   # Handle NA values: remove them before naming or assign them a specific color
#   unique_groups <- unique_groups[!is.na(unique_groups)]
#   unique_families <- unique_families[!is.na(unique_families)]
#
#   # Generate distinct colors for each unique element
#   group_colors <- rainbow(length(unique_groups))
#   names(group_colors) <- unique_groups
#   family_colors <- rainbow(length(unique_families))
#   names(family_colors) <- unique_families
#
#   # Create the column annotation
#   column_ha <- ComplexHeatmap::HeatmapAnnotation(
#     group = KinaseFamily[o, "kinase_group"],
#     family = KinaseFamily[o, "kinase_family"],
#     col = list(
#       group = group_colors,
#       family = family_colors
#     )
#   )
#
#   # Check if sites and o are not empty before drawing the heatmap
#   if (length(sites) > 0 && length(o) > 0) {
#     ht <- ComplexHeatmap::Heatmap(
#       phosScoringMatrices$combinedScoreMatrix[sites, o],
#       top_annotation = column_ha,
#       name = paste("Top", 3, "phosphosite(s)\nfor each kinase")
#     )
#
#     ht <-  ComplexHeatmap::draw(ht)  # Return the Heatmap object
#
#   } else {
#     message("No valid sites or kinases to display in the heatmap for comparison ", comparison)
#     ht <-  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(matrix(0, nrow = 1, ncol = 1), name = "No Data"))
#
#   }
# }
#
# renderHeatmap_DOWN <- function(kinsub_results, comparison) {
#
#   # Read the RDS file reactively
#   phosScoringMatrices <- kinsub_results$kinase_score_mat_list_DOWN[[comparison]]
#
#   # If phosScoringMatrices is NULL, return a placeholder Heatmap
#   if (is.null(phosScoringMatrices) || is.null(phosScoringMatrices$combinedScoreMatrix)) {
#     message("combinedScoreMatrix not found in the data for comparison ", comparison)
#     return(ComplexHeatmap::Heatmap(matrix(0, nrow = 1, ncol = 1), name = "No Data"))
#   }
#
#   # Extract top sites
#   sites <- c()
#   for (i in seq_len(ncol(phosScoringMatrices$combinedScoreMatrix))) {
#     sites <- union(sites, names(
#       sort(phosScoringMatrices$combinedScoreMatrix[, i], decreasing = TRUE)[seq_len(3)]
#     ))
#   }
#
#   o <- intersect(colnames(phosScoringMatrices$combinedScoreMatrix),
#                  rownames(KinaseFamily))
#
#   # Get unique values for group and family
#   unique_groups <- unique(KinaseFamily[o, "kinase_group"])
#   unique_families <- unique(KinaseFamily[o, "kinase_family"])
#
#   # Handle NA values: remove them before naming or assign them a specific color
#   unique_groups <- unique_groups[!is.na(unique_groups)]
#   unique_families <- unique_families[!is.na(unique_families)]
#
#   # Generate distinct colors for each unique element
#   group_colors <- rainbow(length(unique_groups))
#   names(group_colors) <- unique_groups
#   family_colors <- rainbow(length(unique_families))
#   names(family_colors) <- unique_families
#
#   # Create the column annotation
#   column_ha <- ComplexHeatmap::HeatmapAnnotation(
#     group = KinaseFamily[o, "kinase_group"],
#     family = KinaseFamily[o, "kinase_family"],
#     col = list(
#       group = group_colors,
#       family = family_colors
#     )
#   )
#
#   # Check if sites and o are not empty before drawing the heatmap
#   if (length(sites) > 0 && length(o) > 0) {
#     ht <- ComplexHeatmap::Heatmap(
#       phosScoringMatrices$combinedScoreMatrix[sites, o],
#       top_annotation = column_ha,
#       name = paste("Top", 3, "phosphosite(s)\nfor each kinase")
#     )
#
#     ht <-  ComplexHeatmap::draw(ht)  # Return the Heatmap object
#
#   } else {
#     message("No valid sites or kinases to display in the heatmap for comparison ", comparison)
#     ht <-  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(matrix(0, nrow = 1, ncol = 1), name = "No Data"))
#
#   }
# }

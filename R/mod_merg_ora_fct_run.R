#' Run over-representation analysis (ORA) for all comparisons
#'
#' @description
#' Performs over-representation analysis (ORA) on differentially regulated
#' features across all statistical comparisons contained in a \code{stats_data}
#' object. For each comparison, separate enrichment analyses are run for
#' up-regulated and down-regulated features using Gene Ontology (GO),
#' Reactome pathways, and KEGG pathways.
#'
#' The function supports both global proteomics and phosphoproteomics analyses.
#' Gene identifiers are automatically handled according to the selected mode,
#' and enrichment failures or empty gene sets are safely skipped.
#'
#' @param stats_data A list produced by the statistical analysis module.
#'   Must contain at least:
#'   \itemize{
#'     \item \code{results}: a named list of differential analysis result tables,
#'           one per comparison.
#'     \item \code{annotation}: a data.frame containing feature annotation,
#'           including gene symbols.
#'   }
#' @param AnnotDb An \code{AnnotationDb} object (e.g. \code{org.Hs.eg.db})
#'   used for gene identifier mapping and GO enrichment.
#' @param taxo_reactome Character scalar specifying the Reactome organism code
#'   (e.g. \code{"human"}, \code{"mouse"}).
#' @param taxo_kegg Character scalar specifying the KEGG organism code
#'   (e.g. \code{"hsa"}, \code{"mmu"}).
#' @param mode Character scalar indicating the analysis mode. Either
#'   \code{"global"} (protein-level analysis) or \code{"phospho"}
#'   (phosphosite-level analysis).
#'
#' @return
#' A named list with one entry per comparison. Each comparison contains a list
#' of enrichment results with the following elements:
#' \itemize{
#'   \item \code{GO_up}, \code{GO_down}: GO enrichment results for up- and
#'         down-regulated genes.
#'   \item \code{React_up}, \code{React_down}: Reactome pathway enrichment results.
#'   \item \code{KEGG_up}, \code{KEGG_down}: KEGG pathway enrichment results.
#' }
#' Individual elements are \pkg{clusterProfiler} or \pkg{ReactomePA} result
#' objects, or \code{NULL} if no enrichment could be computed.
#'
#' @details
#' Differentially regulated features are identified using the \code{sig.FDR}
#' column, where values of \code{1} indicate up-regulation and \code{-1}
#' indicate down-regulation. Gene symbols are cleaned to retain the first
#' identifier when multiple symbols are present.
#'
#' Enrichment analyses are executed using safe wrappers to prevent errors
#' caused by empty gene sets or mapping failures.
#'
#' @importFrom stats cor na.omit quantile rnorm sd setNames
#' @importFrom AnnotationDbi mapIds
#' @importFrom clusterProfiler enrichGO enrichKEGG
#' @importFrom ReactomePA enrichPathway
#'
#' @export
run_all_ORA <- function(
    stats_data,
    AnnotDb,
    taxo_reactome,
    taxo_kegg,
    mode = c("global", "phospho")
) {

  mode <- match.arg(mode)

  id_col <- "gene"

  message("ORA mode = ", mode)
  message("Using ID column: ", id_col)

  comparisons <- names(stats_data$results)[grepl("_vs_", names(stats_data$results))]

  # ---------- safe wrapper ----------
  safe_enrich <- function(fun, entrezid_vec, ...) {
    if (is.null(entrezid_vec) || length(entrezid_vec) == 0) {
      message("Skipped: no ENTREZID for ", deparse(substitute(fun)))
      return(NULL)
    }

    suppressWarnings(
      tryCatch(
        fun(gene = entrezid_vec, ...),
        error = function(e) {
          message("Enrichment failed: ", deparse(substitute(fun)), " (", e$message, ")")
          NULL
        }
      )
    )
  }

  run_enrichment <- function(comp_name) {

    res <- stats_data$results[[comp_name]] %>%
      tibble::rownames_to_column("uniprot_id") %>%
      dplyr::left_join(stats_data$annotation, by = "uniprot_id")

    # ----- Identify up/down -----
    up_raw   <- unique(res[[id_col]][res$sig.FDR == 1])
    down_raw <- unique(res[[id_col]][res$sig.FDR == -1])

    clean <- function(x) {
      if (length(x) == 0 || all(is.na(x))) return(character(0))
      unique(sapply(strsplit(x, ";"), `[`, 1))
    }

    up_clean <- clean(up_raw)
    down_clean <- clean(down_raw)

    # ----- Safe mapping -----
    safe_map <- function(gene_vec) {
      if (length(gene_vec) == 0) return(NULL)
      AnnotationDbi::mapIds(
        AnnotDb,
        keys = gene_vec,
        keytype = "SYMBOL",
        column = "ENTREZID"
      )
    }

    clean_ids <- function(x) {
      x <- as.character(x)
      x <- x[!is.na(x)]
      unique(x)
    }

    # ----- SYMBOL → ENTREZ -----
    up_entrez   <- clean_ids(safe_map(up_clean))
    down_entrez <- clean_ids(safe_map(down_clean))

    # ----- Enrichments -----
    list(
      GO_up = safe_enrich(clusterProfiler::enrichGO,
                          entrezid_vec = up_entrez,
                          OrgDb = AnnotDb,
                          keyType = "ENTREZID",
                          ont = "ALL",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          readable = TRUE),

      GO_down = safe_enrich(clusterProfiler::enrichGO,
                            entrezid_vec = down_entrez,
                            OrgDb = AnnotDb,
                            keyType = "ENTREZID",
                            ont = "ALL",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            readable = TRUE),

      React_up = safe_enrich(ReactomePA::enrichPathway,
                             entrezid_vec = up_entrez,
                             organism = taxo_reactome,
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05,
                             readable = TRUE),

      React_down = safe_enrich(ReactomePA::enrichPathway,
                               entrezid_vec = down_entrez,
                               organism = taxo_reactome,
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05,
                               readable = TRUE),

      KEGG_up = safe_enrich(clusterProfiler::enrichKEGG,
                            entrezid_vec = up_entrez,
                            organism = taxo_kegg,
                            pvalueCutoff = 0.05),

      KEGG_down = safe_enrich(clusterProfiler::enrichKEGG,
                              entrezid_vec = down_entrez,
                              organism = taxo_kegg,
                              pvalueCutoff = 0.05)
    )
  }

  out <- setNames(lapply(comparisons, run_enrichment), comparisons)
  out
}


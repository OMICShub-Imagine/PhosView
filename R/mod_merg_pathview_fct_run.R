#' Perform KEGG enrichment and prepare data for Pathview visualization
#'
#' Integrates differential phosphosite data and inferred kinase activities to
#' generate KEGG pathway enrichment results and structured input suitable for
#' Pathview visualization.
#'
#' @param stat_data Differential analysis results from \code{phos$result()}.
#'   Must contain \code{$results} and \code{$annotation}.
#' @param kinsub_data Output from \code{kinsub_results()} containing inferred
#'   kinase activity scores.
#' @param AnnotDb An \code{AnnotationDbi} organism database
#'   (e.g. \code{org.Hs.eg.db}).
#' @param taxo_kegg Character string specifying KEGG taxonomy
#'   (e.g. \code{"hsa"}, \code{"mmu"}).
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Builds a unified table of regulated substrates and kinases
#'   \item Performs KEGG over-representation analysis per comparison
#'   \item Maps Entrez IDs back to gene symbols for interpretability
#' }
#'
#' @return A list with two elements:
#' \describe{
#'   \item{final_input}{Combined substrate and kinase input table}
#'   \item{enrichment}{KEGG enrichment results per comparison}
#' }
#'
#' @seealso \code{\link[clusterProfiler]{enrichKEGG}}
#'
#' @export

compute_pathview <- function(
    stat_data,      # phos$result()$stat_data
    kinsub_data,    # kinsub_results()
    AnnotDb,        # org.Hs.eg.db (or other)
    taxo_kegg       # "hsa" / "mmu" / etc
) {

  #==============================================================
  # 1. Build unified table
  #==============================================================

  build_input_table <- function(stat_data, kinsub_data) {

    comparisons <- names(stat_data$results)

    tables <- lapply(comparisons, function(comp) {

      res <- stat_data$results[[comp]] %>%
        tibble::rownames_to_column("phosphosite") %>%
        dplyr::left_join(stat_data$annotation, by = "phosphosite") %>%
        dplyr::filter(sig.FDR %in% c(-1, 1)) %>%
        dplyr::select(phosphosite, gene, logFC, sig.FDR) %>%
        dplyr::mutate(
          Category = dplyr::if_else(sig.FDR == 1, "Substrate_Up", "Substrate_Down"),
          Comparison = comp,
          phosphosite = as.character(phosphosite),
          gene = as.character(gene),
          logFC = as.numeric(logFC),
          sig.FDR = as.numeric(sig.FDR)
        )

      kin_up   <- kinsub_data$kinase_score_mat_list_UP[[comp]]$combinedScoreMatrix
      kin_down <- kinsub_data$kinase_score_mat_list_DOWN[[comp]]$combinedScoreMatrix

      kin_tbls <- list()

      if (!is.null(kin_up) && ncol(kin_up) > 0) {
        kin_tbls[[length(kin_tbls)+1]] <- tibble::tibble(
          phosphosite = NA_character_,
          gene = colnames(kin_up),
          logFC = NA_real_,
          Category = "Kinase_Up",
          Comparison = comp
        )
      }

      if (!is.null(kin_down) && ncol(kin_down) > 0) {
        kin_tbls[[length(kin_tbls)+1]] <- tibble::tibble(
          phosphosite = NA_character_,
          gene = colnames(kin_down),
          logFC = NA_real_,
          Category = "Kinase_Down",
          Comparison = comp
        )
      }

      kinases <- if (length(kin_tbls) > 0) dplyr::bind_rows(kin_tbls) else NULL

      if (!is.null(kinases)) {

        kinases <- kinases %>%
          dplyr::mutate(
            Category = dplyr::case_when(
              grepl("Kinase_Up", Category) &
                gene %in% res$gene[res$Category == "Substrate_Up"]
              ~ "Kinase_Up_Substrate_Up",

              grepl("Kinase_Down", Category) &
                gene %in% res$gene[res$Category == "Substrate_Down"]
              ~ "Kinase_Down_Substrate_Down",

              TRUE ~ Category
            )
          )

        dplyr::bind_rows(res, kinases)
      } else {
        res
      }
    })

    dplyr::bind_rows(tables)
  }

  final_input <- build_input_table(stat_data, kinsub_data)


  if (is.null(final_input) || nrow(final_input) == 0 ||
      all(is.na(final_input$gene))) {

    message("No phosphosite/kinase signal → skipping Pathview pipeline")

    return(list(
      final_input = tibble::tibble(
        phosphosite = character(),
        gene = character(),
        logFC = numeric(),
        sig.FDR = numeric(),
        Category = character(),
        Comparison = character()
      ),
      enrichment = list(),
      status = "empty"
    ))
  }

  #==============================================================
  # 2. KEGG enrichment
  #==============================================================

  comparisons <- unique(final_input$Comparison)

  run_enrichment <- function(comp) {

    genes <- final_input %>%
      dplyr::filter(Comparison == comp) %>%
      dplyr::pull(gene) %>%
      unique() %>%
      stats::na.omit()

    genes <- genes[!is.na(genes) & genes != ""]

    if (length(genes) == 0) return(NULL)

    entrez <- AnnotationDbi::mapIds(
      AnnotDb,
      keys = genes,
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )

    entrez <- entrez[!is.na(entrez)]

    if (length(entrez) == 0) return(NULL)

    KEGG <- tryCatch(
      clusterProfiler::enrichKEGG(
        entrez,
        organism = taxo_kegg,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 500
      ),
      error = function(e) NULL
    )

    if (is.null(KEGG) || nrow(KEGG@result) == 0) return(NULL)

    # =========================
    # keep your function (UNCHANGED)
    # =========================

    map_entrez_to_gene <- function(entrez_str, gene_map) {
      ids <- unlist(strsplit(entrez_str, "/"))
      symbols <- names(gene_map)[match(ids, gene_map)]
      symbols <- symbols[!is.na(symbols)]
      paste(symbols, collapse = "/")
    }

    KEGG@result$geneName <- sapply(
      KEGG@result$geneID,
      map_entrez_to_gene,
      gene_map = entrez
    )

    sig <- KEGG@result %>% dplyr::filter(p.adjust < 0.05)
    if (nrow(sig) == 0)
      sig <- KEGG@result %>% dplyr::filter(pvalue < 0.05)
    if (nrow(sig) == 0) sig <- NULL

    list(KEGG = KEGG, sig = sig)
  }

  enrichment_results <- lapply(comparisons, run_enrichment)
  names(enrichment_results) <- comparisons

  #==============================================================
  # FINAL OUTPUT VALIDATION (SAFE CHECKS ONLY)
  #==============================================================

  validate_pathview_output <- function(res) {

    message("=== PATHVIEW OUTPUT VALIDATION ===")

    # ---- check structure ----
    if (is.null(res)) stop("compute_pathview returned NULL")

    if (!is.list(res)) stop("Output is not a list")

    expected_names <- c("final_input", "enrichment")

    missing <- setdiff(expected_names, names(res))
    if (length(missing) > 0) {
      stop("Missing top-level elements: ", paste(missing, collapse = ", "))
    }

    message("Top-level structure OK")

    # ---- final_input checks ----
    fi <- res$final_input

    if (is.null(fi)) stop("final_input is NULL")
    if (!is.data.frame(fi)) stop("final_input is not a data.frame")

    message("final_input is data.frame")

    required_cols <- c("gene", "Category", "Comparison")

    missing_cols <- setdiff(required_cols, names(fi))
    if (length(missing_cols) > 0) {
      stop("final_input missing columns: ", paste(missing_cols, collapse = ", "))
    }

    message("final_input columns OK")

    if (nrow(fi) == 0) stop("final_input has 0 rows")

    message("final_input rows: ", nrow(fi))

    # ---- enrichment checks ----
    en <- res$enrichment

    if (is.null(en)) warning("enrichment is NULL")

    if (is.list(en)) {
      message("enrichment is a list with ", length(en), " elements")

      comp_names <- names(en)
      message("comparisons: ", paste(comp_names, collapse = ", "))
    }

    message("=== VALIDATION COMPLETE ===")

    invisible(TRUE)
  }

  res <- list(
    final_input = final_input,
    enrichment  = enrichment_results
  )

  validate_pathview_output(res)

  return(res)
}


#' Generate a Pathview KEGG pathway PNG file
#'
#' Runs \code{pathview()} for a given KEGG pathway and returns the file path
#' to the generated PNG image.
#'
#' @param geneList Named numeric vector of gene values (Entrez IDs as names).
#' @param pathway_id Character KEGG pathway identifier (e.g. \code{"hsa04010"}).
#' @param taxo_kegg Character KEGG organism code (\code{"hsa"} or \code{"mmu"}).
#' @param out_dir Output directory where the PNG file will be written.
#' @param suffix Character suffix appended to the output file name.
#'
#' @return Character string giving the full path to the generated PNG file.
#'
#' @seealso \code{\link[pathview]{pathview}}
#'
#' @export

generate_pathview_png <- function(
  geneList,
  pathway_id,
  taxo_kegg,
  out_dir,
  suffix
) {

  stopifnot(taxo_kegg %in% c("hsa", "mmu"))
  stopifnot(!is.na(pathway_id))
  stopifnot(length(pathway_id) == 1)
  stopifnot(pathway_id != "")

  cat("\n=== GENERATE PATHVIEW ===\n")
  cat("pathway_id =", pathway_id, "\n")
  cat("taxo =", taxo_kegg, "\n")
  cat("out_dir =", out_dir, "\n")
  cat("gene count =", length(geneList), "\n")


  withr::with_dir(out_dir,
                  pathview::pathview(
                  gene.data   = geneList,
                  pathway.id  = pathway_id,
                  species     = taxo_kegg,
                  gene.idtype = "entrez",
                  limit      = list(gene = c(-max(abs(geneList)), max(abs(geneList))), cpd = 1),
                  kegg.native = TRUE,
                  out.suffix  = suffix,
                  kegg.dir    = out_dir,
                  low  = list(gene = "green"),
                  mid  = list(gene = "gray"),
                  high = list(gene = "red")
  ))

  cat("\nFiles after pathview:\n")
  print(list.files(out_dir, full.names = TRUE))

  files <- list.files(
    out_dir,
    pattern = paste0(pathway_id, ".*", suffix, ".*\\.png$"),
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("Pathview PNG not found in: ", out_dir)
  }

  files[1]
}

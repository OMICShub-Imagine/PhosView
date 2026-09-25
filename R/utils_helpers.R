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

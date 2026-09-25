#' Load and validate uploaded quantitative or metadata files
#'
#' This function loads a tab-delimited file uploaded by the user and performs
#' basic validation depending on the selected data type. It supports:
#' \itemize{
#'   \item Quantitative phosphoproteomics output files (e.g. Spectronaut)
#'   \item Experimental design / metadata files
#' }
#'
#' @param file Character scalar giving the full path to the uploaded file.
#'   The file must exist and be readable.
#'
#' @param data_type Character indicating the type of data to load.
#'   Must be one of:
#'   \describe{
#'     \item{\code{"quant"}}{Quantitative phosphoproteomics data file.
#'       The file is expected to be tab-delimited and to contain a predefined
#'       set of protein- and PTM-level annotation columns.}
#'     \item{\code{"meta"}}{Experimental metadata file describing samples,
#'       conditions, and groups.}
#'   }
#'
#' @details
#' When \code{data_type = "quant"}, the function:
#' \itemize{
#'   \item Reads the file using \code{\link[utils]{read.table}}
#'   \item Disables automatic name repair to preserve original column names
#'   \item Verifies the presence of mandatory protein and PTM annotation columns
#' }
#'
#' Required columns for quantitative data include (non-exhaustive):
#' \itemize{
#'   \item \code{PG.ProteinGroups}
#'   \item \code{PG.Genes}
#'   \item \code{PTM.CollapseKey}
#'   \item \code{PTM.ModificationTitle}
#'   \item \code{PTM.SiteAA}
#'   \item \code{PTM.SiteLocation}
#' }
#'
#' When \code{data_type = "meta"}, the function:
#' \itemize{
#'   \item Reads a tab-delimited experimental design file
#'   \item Checks that at least three columns are present
#'   \item Returns the metadata unchanged for downstream validation
#' }
#'
#' @return
#' A \code{data.frame} containing the validated quantitative data or metadata.
#'
#' @section Errors:
#' The function stops with an error if:
#' \itemize{
#'   \item The file path is NULL, empty, or does not exist
#'   \item Required columns are missing from the quantitative file
#'   \item The metadata file has fewer than three columns
#' }
#'
#' @importFrom utils read.table
#'
#' @examples
#' \dontrun{
#' ## Load quantitative data
#' quant <- mod_phos_first_class_browse_fct_load(
#'   file = "phospho_quant.tsv",
#'   data_type = "quant"
#' )
#'
#' ## Load experimental metadata
#' meta <- mod_phos_first_class_browse_fct_load(
#'   file = "exp_design.tsv",
#'   data_type = "meta"
#' )
#' }
#'
#' @export
mod_phos_first_class_browse_fct_load <- function(file, data_type = c("quant", "meta")) {
  data_type <- match.arg(data_type)

  if (is.null(file) || !nzchar(file) || !file.exists(file)) {
    stop("`file` is NULL or does not exist: ", as.character(file))
  }

  if (data_type == "quant") {
    message("Reading quantitative data from: ", file)

    input <- read.table(
      file = file,
      dec = ",",
      sep = "\t",
      header = TRUE,
      quote = "",
      check.names = FALSE
    )

    required_cols <- c(
      "PG.ProteinGroups",
      "PG.Genes",
      "PG.Organisms",
      "PG.ProteinDescriptions",
      "PG.ProteinNames",
      "PG.NrOfStrippedSequencesIdentified (Experiment-wide)",
      "PG.NrOfModifiedSequencesIdentified (Experiment-wide)",
      "PG.NrOfPrecursorsIdentified (Experiment-wide)",
      "PG.Completeness",
      "PG.CV",
      "PTM.ProteinId",
      "PTM.CollapseKey",
      "PTM.Multiplicity",
      "PTM.ModificationTitle",
      "PTM.SiteAA",
      "PTM.SiteLocation",
      "PTM.FlankingRegion"

    )

    missing_cols <- setdiff(required_cols, colnames(input))
    if (length(missing_cols) > 0) {
      stop(
        "The following required column(s) are missing in the quantitative file: ",
        paste(missing_cols, collapse = ", ")
      )
    }

    message("Quantitative data loaded and validated.")
    return(as.data.frame(input, stringsAsFactors = FALSE))
  }

  if (data_type == "meta") {
    message("Reading experimental metadata from: ", file)

    meta <- read.table(
      file = file,
      sep = "\t",
      header = TRUE
    ) %>% as.data.frame()


    if (ncol(meta) < 3) {
      stop("Experimental metadata file must have at least 3 columns (Sample_ID, Simplify_ID, group).")
    }

    message("Experimental metadata loaded (", nrow(meta), " rows, ", ncol(meta), " cols).")
    return(meta)
  }
}

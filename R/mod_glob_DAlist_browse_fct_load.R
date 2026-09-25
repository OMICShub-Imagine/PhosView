#' Load and validate an uploaded quantitative or metadata file
#'
#' @description
#' Reads a tab-delimited file containing either quantitative protein group data
#' or experimental metadata, performs basic validation checks, and returns the
#' content as a \code{data.frame}.
#'
#' @details
#' The behavior depends on the value of \code{data_type}:
#' \itemize{
#'   \item \strong{"quant"}: The file is expected to be a tab-separated quantitative
#'   data table (e.g. DIA-NN or Spectronaut protein group output with specific columns (cf vignette)). A predefined set of required
#'   columns is checked, and the function stops with an error if any are missing.
#'
#'   \item \strong{"meta"}: The file is expected to be a tab-separated experimental
#'   metadata table. It must contain at least three columns (typically corresponding
#'   to Sample_ID, Simplified_ID, and experimental group).
#' }
#'
#' The function stops with an informative error message if the file does not exist,
#' is empty, or does not meet the expected structure.
#'
#' @param file Character scalar. Full path to the uploaded tab-delimited file.
#'   The file must exist and be readable.
#' @param data_type Character scalar. Type of data to load.
#'   Must be one of \code{"quant"} (quantitative data) or \code{"meta"}
#'   (experimental metadata).
#'
#' @return
#' A \code{data.frame} containing the validated input data.
#' For quantitative data, all required columns are guaranteed to be present.
#'
#' @importFrom utils data head read.table str
#'
#' @examples
#' \dontrun{
#' ## Load quantitative data
#' quant_df <- mod_glob_DAlist_browser_fct_load(
#'   file = "proteinGroups.txt",
#'   data_type = "quant"
#' )
#'
#' ## Load experimental metadata
#' meta_df <- mod_glob_DAlist_browser_fct_load(
#'   file = "metadata.tsv",
#'   data_type = "meta"
#' )
#' }
#'
#' @export
mod_glob_DAlist_browser_fct_load <- function(file, data_type = c("quant", "meta")) {
  data_type <- match.arg(data_type)

  if (is.null(file) || !nzchar(file) || !file.exists(file)) {
    stop("`file` is NULL or does not exist: ", as.character(file))
  }

  if (data_type == "quant") {
    message("Reading quantitative data from: ", file)

    probe <- read.table(
      file = file,
      sep = "\t",
      header = TRUE,
      quote = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    # ---- DIA-NN ----
    if ("Protein.Group" %in% colnames(probe)) {

      message("Detected DIA-NN output")

      input <- read.table(
        file = file,
        dec = ".",
        sep = "\t",
        header = TRUE,
        quote = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      required_cols <- c(
        "Protein.Group",
        "Protein.Names",
        "Genes",
        "First.Protein.Description",
        "N.Sequences",
        "N.Proteotypic.Sequences"
      )

      # ---- Spectronaut ----
    } else if ("PG.ProteinGroups" %in% colnames(probe)) {

      message("Detected Spectronaut output")

      input <- read.table(
        file = file,
        dec = ",",
        sep = "\t",
        header = TRUE,
        quote = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
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
        "PG.MolecularWeight"
      )

      # ---- UNKNOWN ----
    } else {
      stop("Unknown quantitative file format: could not detect required columns.")
    }

    # ---- COMMON VALIDATION ----
    missing_cols <- setdiff(required_cols, colnames(input))

    if (length(missing_cols) > 0) {
      stop(
        "Missing required column(s): ",
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
      header = TRUE,
      stringsAsFactors = FALSE
    )

    meta <- as.data.frame(meta, stringsAsFactors = FALSE)

    if (ncol(meta) < 3) {
      stop("Experimental metadata file must have at least 3 columns (Sample_ID, Simplify_ID, group).")
    }

    message("Experimental metadata loaded (", nrow(meta), " rows, ", ncol(meta), " cols).")
    return(meta)
  }
}

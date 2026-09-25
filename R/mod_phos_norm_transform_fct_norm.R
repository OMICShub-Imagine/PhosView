#' Normalize phosphosite DAList data
#'
#' This function normalizes phosphosite intensity data contained in a \code{DAList} object.
#' Several normalization methods are available, all followed by a log2 transformation. After
#' log2 transformation, a +10 shift is applied to avoid negative values (useful for Perseus-like
#' downstream processing).
#'
#' @param DAList_obj A \code{DAList} object from \pkg{proteoDA} containing phosphosite intensity
#'   data in the \code{$data} slot.
#' @param norm Character. The normalization method to apply. Choices include:
#'   \itemize{
#'     \item "no normalisation and log2 transformation"
#'     \item "median normalisation and log2 transformation"
#'     \item "mean normalisation and log2 transformation"
#'     \item "vsn (variance stabilizing normalization) normalisation and log2 transformation"
#'     \item "quantile normalisation and log2 transformation"
#'     \item "cycloess (cyclic loess) normalisation and log2 transformation"
#'     \item "rlr (robust linear regression) normalisation and log2 transformation"
#'     \item "gi (global intensity) normalisation and log2 transformation"
#'   }
#'
#' @details
#' The function internally maps the user-friendly selection to the appropriate method supported
#' by \pkg{proteoDA}'s \code{normalize_data} function. After normalization, all values are log2
#' transformed and shifted by +10 to prevent negative values, which can cause issues in some
#' downstream analyses such as Perseus-like imputation.
#'
#' @return A normalized \code{DAList} object with transformed data in the \code{$data} slot.
#'
#' @importFrom proteoDA normalize_data
#' @export
#'
#' @examples
#' \dontrun{
#' # Median normalization and log2 transformation
#' DAList_norm <- mod_phos_norm_transform_fct_norm(
#'   DAList_obj = my_DAList,
#'   norm = "median normalisation and log2 transformation"
#' )
#' }
mod_phos_norm_transform_fct_norm <- function(DAList_obj, norm) {

  # Basic checks
  if (is.null(DAList_obj)) stop("DAlist_object is NULL.")
  if (missing(norm) || length(norm) == 0) {
    warning("No normalization provided: returning original DAlist.")
    return(DAList_obj)
  }

  # Convert text selection to correct argument function
  norm_arg <- switch(norm,
                     "no normalization but log2 transformation" = "log2",
                     "median normalization followed by log2 transformation" = "median",
                     "mean normalization followed by log2 transformation" = "mean",
                     "vsn (variance stabilizing normalization) normalization followed by log2 transformation" = "vsn",
                     "quantile normalization followed by log2 transformation" = "quantile",
                     "cycloess (cyclic loess) normalization followed by log2 transformation" = "cycloess",
                     "rlr (robust linear regression) normalization followed by log2 transformation" = "rlr",
                     "gi (global intensity) normalization followed by log2 transformation" = "gi",
                     stop("Unknown normalization method selection"))


  # Call proteoDA::filter_samples correctly
  norm_data <- proteoDA::normalize_data(
    DAList = DAList_obj,
    norm_method = norm_arg
  )

  ## Add +10 in all log2 transform value (to avoid possible negative value with Perseus-like Imputation)
  norm_data$data <- norm_data$data + 10

  return(norm_data)
}

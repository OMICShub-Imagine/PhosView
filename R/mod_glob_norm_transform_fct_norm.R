#' Normalize and log-transform quantitative data in a DAList object
#'
#' @description
#' Applies normalization and log2 transformation to the quantitative data stored
#' in a \code{proteoDA::DAList} object using one of several predefined
#' normalization methods. The function acts as a wrapper around
#' \code{proteoDA::normalize_data}.
#'
#' @details
#' The normalization method is selected via a human-readable character string,
#' typically originating from a user interface. This selection is internally
#' mapped to the corresponding \code{norm_method} argument expected by
#' \code{proteoDA::normalize_data}.
#'
#' After normalization and log2 transformation, a constant value of \code{10} is
#' added to all quantitative values to avoid negative values, which can be
#' problematic for downstream steps such as Perseus-like imputation.
#'
#' If no normalization method is provided, the input \code{DAList} object is
#' returned unchanged with a warning.
#'
#' @param DAList_obj A \code{proteoDA::DAList} object containing quantitative
#'   data to be normalized.
#' @param norm Character scalar specifying the normalization method. Supported
#'   values include:
#'   \itemize{
#'     \item \code{"no normalisation and log2 transformation"}
#'     \item \code{"median normalisation and log2 transformation"}
#'     \item \code{"mean normalisation and log2 transformation"}
#'     \item \code{"vsn (variance stabilizing normalization) normalisation and log2 transformation"}
#'     \item \code{"quantile normalisation and log2 transformation"}
#'     \item \code{"cycloess (cyclic loess) normalisation and log2 transformation"}
#'     \item \code{"rlr (robust linear regression) normalisation and log2 transformation"}
#'     \item \code{"gi (global intensity) normalisation and log2 transformation"}
#'   }
#'
#' @return
#' A \code{proteoDA::DAList} object with normalized and log2-transformed data.
#'
#' @seealso
#' \code{\link[proteoDA]{normalize_data}}
#'
#' @examples
#' \dontrun{
#' ## Apply median normalization
#' da_norm <- mod_glob_norm_transform_fct_norm(
#'   DAList_obj = da_list,
#'   norm = "median normalisation and log2 transformation"
#' )
#' }
#'
#' @export
mod_glob_norm_transform_fct_norm <- function(DAList_obj, norm) {

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

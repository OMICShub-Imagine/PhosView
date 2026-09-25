#' Filter samples from a DAList object using a logical condition
#'
#' This internal helper filters samples from a \code{proteoDA::DAList} object
#' based on a logical condition evaluated on the metadata. Samples for which
#' the condition evaluates to \code{TRUE} are retained, while others are removed.
#'
#' The filtering condition is evaluated within the metadata context and the
#' calling environment, allowing flexible, expression-based filtering.
#'
#' @param DAList A \code{proteoDA::DAList} object containing quantitative data
#'   and sample metadata.
#'
#' @param condition A logical expression evaluated in the context of
#'   \code{DAList$metadata}. Samples satisfying the condition are kept.
#'   For example: \code{!Simplify_ID \%in\% outliers}.
#'
#' @details
#' The function updates both the \code{metadata} and \code{data} slots of the
#' \code{DAList} object to ensure consistency. Informative messages are printed
#' to report removed samples.
#'
#' This function is adapted from \code{proteoDA::filter_samples} with extended
#' support for evaluating conditions in the caller environment.
#'
#' @return A \code{proteoDA::DAList} object with filtered samples.
#'
#' @export

filter_samples2 <- function(DAList, condition) {

  if (!(class(DAList) %in% c("DAList"))) {
    cli::cli_abort("{.arg DAList} must be a DAList object")
  }

  if (is.null(DAList$metadata)) {
    cli::cli_abort("{.arg DAList} does not contain metadata for filtering samples")
  }

  # get input metadata
  in_meta <- DAList$metadata

  # Capture the expression
  condition_call <- substitute(condition)

  # Evaluate condition in both metadata and caller environment (modified part from original proteoDA::filter_samples function)
  caller_env <- parent.frame()
  in_meta <- within(in_meta, {
    sample_to_be_kept <- eval(condition_call, envir = c(as.list(in_meta), as.list(environment())), enclos = caller_env)
  })

  # Subset as before
  meta_removed <- subset(in_meta, subset = !sample_to_be_kept, select = c(-sample_to_be_kept))
  meta_kept <- subset(in_meta, subset = sample_to_be_kept, select = c(-sample_to_be_kept))

  if (nrow(meta_removed) + nrow(meta_kept) != nrow(in_meta)) {
    cli::cli_abort(c("Issue when filtering samples",
                     "!" = "Rows kept + rows removed != rows input",
                     "i" = "Is there an {.val NA} in the column you are filtering on?"))
  }

  DAList$metadata <- meta_kept
  DAList$data <- DAList$data[, rownames(DAList$metadata)]
  proteoDA:::validate_DAList(DAList) ### modify adding "proteoDA::" if not it will not be found

  cli::cli_inform("Removed {.val {nrow(meta_removed)}} of the {.val {nrow(in_meta)}} {cli::qty(nrow(in_meta))} sample{?s} in DAList")
  cli::cli_inform("{cli::qty(nrow(meta_removed))} Sample{?s} removed:")
  message(paste0(utils::capture.output(meta_removed), collapse = "\n"))

  DAList
}


#' Filter outlier samples from a phosphoproteomics DAList object
#'
#' This helper function removes user-specified outlier samples from a
#' \code{proteoDA::DAList} object based on their \code{Simplify_ID}.
#' It is designed as a high-level wrapper for interactive or pipeline-based
#' phosphoproteomics workflows.
#'
#' @param DAList_obj A \code{proteoDA::DAList} object containing phosphosite-level
#'   quantitative data and metadata.
#'
#' @param outliers Character vector of \code{Simplify_ID} values corresponding
#'   to samples that should be removed.
#'
#' @details
#' The function applies sample filtering using \code{filter_samples2()}, keeping
#' only samples whose \code{Simplify_ID} is not listed in \code{outliers}.
#'
#' If no outliers are provided, the input \code{DAList} object is returned
#' unchanged.
#'
#' @return A \code{proteoDA::DAList} object with outlier samples removed.
#'
#' @seealso \code{\link{filter_samples2}}
#'
#' @export

mod_phos_filt_sample_transform_fct_filt_samples <- function(DAList_obj, outliers) {

  # Basic checks
  if (is.null(DAList_obj)) stop("DAlist_object is NULL.")
  if (missing(outliers) || length(outliers) == 0) {
    warning("No outliers provided: returning original DAlist.")
    return(DAList_obj)
  }
  if (!"Simplify_ID" %in% names(DAList_obj$metadata)) {
    stop("Simplify_ID column not found in DAlist metadata.")
  }


  # Call proteoDA::filter_samples correctly
  filt_samples <- filter_samples2(
    DAList = DAList_obj,
    condition = !Simplify_ID %in% outliers
  )

  return(filt_samples)
}

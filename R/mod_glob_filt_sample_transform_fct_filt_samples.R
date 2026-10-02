#' Filter samples from a DAList object using a logical condition
#'
#' @description
#' Removes samples from a \code{proteoDA::DAList} object based on a logical
#' condition evaluated on the metadata. This function is adapted from
#' \code{proteoDA::filter_samples} with extended evaluation to allow variables
#' from the calling environment.
#'
#' @param DAList A valid \code{proteoDA::DAList} object containing quantitative
#'   data and metadata.
#' @param condition A logical expression evaluated in the context of the
#'   metadata. Samples for which the condition evaluates to \code{TRUE} are
#'   kept; others are removed.
#'
#' @return
#' A \code{proteoDA::DAList} object with samples filtered according to
#' \code{condition}.
#'
#' @details
#' The condition is evaluated using non-standard evaluation, allowing direct
#' reference to metadata columns (e.g. \code{Simplify_ID}, \code{group}) as well
#' as objects from the calling environment.
#'
#' Informative messages are printed indicating how many samples were removed
#' and which samples were affected.
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
  validate_DAList_local(DAList)

  cli::cli_inform("Removed {.val {nrow(meta_removed)}} of the {.val {nrow(in_meta)}} {cli::qty(nrow(in_meta))} sample{?s} in DAList")
  cli::cli_inform("{cli::qty(nrow(meta_removed))} Sample{?s} removed:")
  message(paste0(utils::capture.output(meta_removed), collapse = "\n"))

  DAList
}


#' Remove outlier samples from a DAList object using Simplify_ID
#'
#' @description
#' Filters out samples from a \code{proteoDA::DAList} object whose
#' \code{Simplify_ID} values are listed as outliers. This function acts as a
#' wrapper around \code{filter_samples2}.
#'
#' @param DAList_obj A \code{proteoDA::DAList} object containing quantitative
#'   data and metadata.
#' @param outliers Character vector of \code{Simplify_ID} values corresponding
#'   to samples that should be removed.
#'
#' @return
#' A \code{proteoDA::DAList} object with outlier samples removed.
#'
#' @examples
#' \dontrun{
#' da_filt <- mod_glob_filt_sample_transform_fct_filt_samples(
#'   DAList_obj = da_list,
#'   outliers = c("S01", "S07")
#' )
#' }
#'
#' @export
mod_glob_filt_sample_transform_fct_filt_samples <- function(DAList_obj, outliers) {

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

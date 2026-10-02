#' Filter proteins based on proportion of valid values in at least one experimental groups
#'
#' @description
#' Filters proteins from a \code{proteoDA::DAList} object based on the proportion
#' of non-missing quantitative values within at least one experimental group.
#' A protein is retained if it meets the minimum proportion threshold in
#' \strong{any} group.
#'
#' @details
#' For each experimental group defined in the metadata, the function calculates
#' the minimum number of non-missing values required for a protein to pass the
#' filtering criterion. Proteins that satisfy this requirement in at least one
#' group are retained.
#'
#' Filtering metadata is recorded in the \code{tags} slot of the returned
#' \code{DAList} object.
#'
#' @param DAList A valid \code{proteoDA::DAList} object containing quantitative
#'   data, annotation, and metadata.
#' @param min_prop Numeric scalar between 0 and 1. Minimum proportion of
#'   non-missing values required within at least one group.
#'   For example, \code{0.5} requires 50\% of samples in a group to have valid
#'   values.
#' @param grouping_column Character scalar. Name of the column in the metadata
#'   slot defining experimental groups.
#'
#' @return
#' A filtered \code{proteoDA::DAList} object containing only proteins that meet
#' the specified proportion threshold in at least one group.
#'
#' @importFrom cli cli_abort cli_inform
#'
#' @examples
#' \dontrun{
#' ## Keep proteins present in at least 50% of samples in one group
#' da_filt <- filter_proteins_by_proportion_onegrp(
#'   DAList = da_list,
#'   min_prop = 0.5,
#'   grouping_column = "group"
#' )
#' }
#'
#' @export
filter_proteins_by_proportion_onegrp <- function(DAList,
                                                 min_prop = NULL,
                                                 grouping_column = "group") {
  validate_DAList_local(DAList)

  if (is.null(DAList$metadata)) {
    cli::cli_abort(c("The {.arg DAList} object must include metadata to filter proteins by group"))
  }

  # Make sure the group column is present in the metadata
  if (!(grouping_column %in% colnames(DAList$metadata))) {
    cli::cli_abort(c("Column {.arg {grouping_column}} not found in metadata slot of {.arg DAList}"))
  }
  # Make sure that min_prop is specified and within the valid range
  if (is.null(min_prop)) {
    cli::cli_abort(c("You must specify a value for {.arg min_prop}",
                     "i" = "A good rule of thumb is to use 0.5 for 50%"))
  }

  if (min_prop <= 0 | min_prop > 1) {
    cli::cli_abort(c("{.arg min_prop} must be from 0 to 1, not {.val {min_prop}}"))
  }

  tmpData <- DAList$data

  # Prep for filtering by getting a vector of group memberships per sample
  # and calculating the threshold for each group
  group_membership <- as.character(DAList$metadata[, grouping_column])

  # Build a results data frame to store whether each protein passes the threshold for each group
  protein_passes_threshold_per_group <- as.data.frame(
    matrix(
      data = NA,
      nrow = nrow(tmpData),
      ncol = length(unique(group_membership)),
      dimnames = list(rownames(tmpData),
                      unique(group_membership))
    )
  )

  # Loop over groups, check if protein passes for that group
  for (group in unique(group_membership)) {
    one_group_data <- tmpData[, group_membership == group]
    # Calculate the minimum number of non-missing values required in this group for a protein to pass
    group_threshold <- floor(ncol(one_group_data) * min_prop)
    # Check if each protein has at least group_threshold non-missing values in this group
    protein_passes_threshold_per_group[, group] <- rowSums(!is.na(one_group_data)) >= group_threshold
  }

  # Keep proteins that pass the threshold in at least one group
  protein_meets_threshold <- apply(X = protein_passes_threshold_per_group, MARGIN = 1, FUN = any)
  # Filter
  kept_proteins <- tmpData[protein_meets_threshold, ]
  removed_proteins <- tmpData[!protein_meets_threshold, ]

  cli::cli_inform("Filtered {.val {nrow(removed_proteins)}} entr{?y/ies} {cli::qty(nrow(removed_proteins))} from the dataset leaving {.val {nrow(kept_proteins)}} entr{?y/ies} {cli::qty(nrow(kept_proteins))} for analysis")

  out <- DAList
  # Update data and annotation
  out$data <- kept_proteins
  out$annotation <- out$annotation[rownames(out$data), , drop = F]
  # Add tags to track filtering
  out$tags$filter_proteins_by_proportion <- c(out$tags$filter_proteins_by_proportion, list(min_prop = min_prop, grouping_column = grouping_column))

  validate_DAList_local(out)

  return(out)
}

#' Filter proteins in a DAList object using predefined proportion thresholds
#'
#' @description
#' Applies protein filtering to a \code{proteoDA::DAList} object based on a
#' user-selected textual threshold describing the minimum proportion of valid
#' values required in at least one experimental group.
#'
#' @details
#' This function serves as a wrapper around
#' \code{\link{filter_proteins_by_proportion_onegrp}}. It converts a
#' human-readable threshold description (e.g. \code{"50% of VV in at least one of
#' your experimental group"}) into a numeric proportion and applies the
#' corresponding filtering.
#'
#' If no threshold is provided, the original \code{DAList} object is returned
#' unchanged.
#'
#' @param DAList_obj A valid \code{proteoDA::DAList} object containing quantitative
#'   data, annotation, and metadata.
#' @param threshold Character scalar describing the filtering threshold.
#'   Must correspond to one of the predefined proportion levels
#'   (50, 70, or 100 percent of valid values in at least one experimental group).
#'
#' @return
#' A \code{proteoDA::DAList} object filtered according to the selected threshold.
#'
#' @examples
#' \dontrun{
#' ## Apply a 70% valid value threshold
#' da_filt <- mod_glob_filt_prot_transform_fct_filt_prot(
#'   DAList_obj = da_list,
#'   threshold = "70% of VV in at least one of your experimental group"
#' )
#' }
#'
#' @export
mod_glob_filt_prot_transform_fct_filt_prot <- function(DAList_obj, threshold) {

  # Basic checks
  if (is.null(DAList_obj)) stop("DAlist_object is NULL.")
  if (missing(threshold) || length(threshold) == 0) {
    warning("No threshold provided, returning original DAlist.")
    return(DAList_obj)
  }

  # Convert text selection to numeric proportion
  min_prop <- switch(threshold,
                     "50% of VV in at least one of your experimental group" = 0.5,
                     "70% of VV in at least one of your experimental group" = 0.7,
                     "100% of VV in at least one of your experimental group" = 1.0,
                     stop("Unknown threshold selection"))

  # Call your helper function
  filt_prot <- filter_proteins_by_proportion_onegrp(
    DAList = DAList_obj,
    min_prop = min_prop,
    grouping_column = "group"
  )

  return(filt_prot)
}

#' Filter phosphosites by valid value proportion in at least one group
#'
#' This internal helper filters phosphosites from a
#' \code{proteoDA::DAList} object by retaining only those that have a minimum
#' proportion of non-missing (valid) quantitative values in at least one
#' experimental group.
#'
#' The filtering is applied independently within each group defined in the
#' metadata. A phosphosite is kept if it satisfies the minimum proportion
#' criterion in at least one group.
#'
#' @param DAList A \code{proteoDA::DAList} object containing phosphosite-level
#'   quantitative data, annotation, and metadata.
#'
#' @param min_prop Numeric value between 0 and 1 specifying the minimum
#'   proportion of valid (non-missing) values required within a group.
#'   For example, \code{0.5} corresponds to 50\%.
#'
#' @param grouping_column Character string specifying the column in
#'   \code{DAList$metadata} that defines experimental groups.
#'   Defaults to \code{"group"}.
#'
#' @details
#' For each experimental group, the function computes the minimum number of
#' non-missing values required for a phosphosite to pass the filter. Phosphosites
#' passing the threshold in at least one group are retained.
#'
#' Filtering information is stored in the \code{tags} slot of the returned
#' \code{DAList} object to ensure traceability.
#'
#' @return A filtered \code{proteoDA::DAList} object containing only phosphosites
#'   meeting the minimum valid value proportion criterion.
#'
#' @export

filter_phosphosite_by_proportion_onegrp <- function(DAList,
                                                 min_prop = NULL,
                                                 grouping_column = "group") {
  proteoDA:::validate_DAList(DAList)

  if (is.null(DAList$metadata)) {
    cli::cli_abort(c("The {.arg DAList} object must include metadata to filter phosphosite by group"))
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
  phosphosite_passes_threshold_per_group <- as.data.frame(
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
    phosphosite_passes_threshold_per_group[, group] <- rowSums(!is.na(one_group_data)) >= group_threshold
  }

  # Keep proteins that pass the threshold in at least one group
  phosphosite_meets_threshold <- apply(X = phosphosite_passes_threshold_per_group, MARGIN = 1, FUN = any)
  # Filter
  kept_phosphosite <- tmpData[phosphosite_meets_threshold, ]
  removed_phosphosite <- tmpData[!phosphosite_meets_threshold, ]

  cli::cli_inform("Filtered {.val {nrow(removed_phosphosite)}} entr{?y/ies} {cli::qty(nrow(removed_phosphosite))} from the dataset leaving {.val {nrow(kept_phosphosite)}} entr{?y/ies} {cli::qty(nrow(kept_phosphosite))} for analysis")

  out <- DAList
  # Update data and annotation
  out$data <- kept_phosphosite
  out$annotation <- out$annotation[rownames(out$data),, drop = F]
  # Add tags to track filtering
  out$tags$filter_phosphosite_by_proportion <- c(out$tags$filter_phosphosite_by_proportion, list(min_prop = min_prop, grouping_column = grouping_column))

  proteoDA:::validate_DAList(out)

  return(out)
}

#' Filter phosphosites based on valid value proportion in at least one group
#'
#' This helper filters phosphosites from a \code{proteoDA::DAList} object by
#' retaining only those with a minimum proportion of valid (non-missing) values
#' in at least one experimental group.
#'
#' @param DAList_obj A \code{proteoDA::DAList} object containing phosphosite-level
#'   quantitative data.
#'
#' @param threshold Character scalar describing the filtering threshold.
#'   Must correspond to one of the predefined proportion levels
#'   (50, 70, or 100 percent of valid values in at least one experimental group).
#'
#' @details
#' Internally, the selected threshold is converted to a numeric proportion
#' (\code{0.5}, \code{0.7}, or \code{1.0}) and applied using
#' \code{filter_phosphosite_by_proportion_onegrp()} with grouping based on the
#' \code{group} column in the metadata.
#'
#' If no threshold is provided, the input \code{DAList} object is returned
#' unchanged.
#'
#' @return A filtered \code{proteoDA::DAList} object containing only phosphosites
#'   meeting the valid value proportion criterion.
#'
#' @seealso \code{\link{filter_phosphosite_by_proportion_onegrp}}
#'
#' @export

mod_phos_filt_phos_transform_fct_filt_phos <- function(DAList_obj, threshold) {

  # Basic checks
  if (is.null(DAList_obj)) stop("DAlist_object is NULL.")
  if (missing(threshold) || length(threshold) == 0) {
    warning("No threshold provided: returning original DAlist.")
    return(DAList_obj)
  }

  # Convert text selection to numeric proportion
  min_prop <- switch(threshold,
                     "50% of VV in at least one of your experimental group" = 0.5,
                     "70% of VV in at least one of your experimental group" = 0.7,
                     "100% of VV in at least one of your experimental group" = 1.0,
                     stop("Unknown threshold selection"))

  # Call your helper function
  filt_phos <- filter_phosphosite_by_proportion_onegrp(
    DAList = DAList_obj,
    min_prop = min_prop,
    grouping_column = "group"
  )

  return(filt_phos)
}

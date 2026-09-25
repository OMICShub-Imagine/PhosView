#' Impute missing values in a DAList object using a left-shifted normal distribution
#'
#' @description
#' Performs missing value imputation on the quantitative data stored in a
#' \code{proteoDA::DAList} object. Missing values are replaced by random draws
#' from a normal distribution that is left-shifted and scaled relative to the
#' observed data distribution of each sample.
#'
#' @details
#' For each sample (column) in the quantitative data matrix:
#' \itemize{
#'   \item The mean, median, and standard deviation of non-missing values are
#'   computed.
#'   \item Missing values are imputed using \code{stats::rnorm()} with:
#'   \itemize{
#'     \item mean = median − (\code{shift} × SD)
#'     \item sd = SD × \code{scale}
#'   }
#' }
#'
#' This approach is commonly used in proteomics to simulate low-abundance values
#' below the detection limit.
#'
#' The function stops with an error if no missing values are present in the data.
#'
#' @param DAlist_obj A \code{proteoDA::DAList} object containing quantitative
#'   data in its \code{data} slot.
#' @param scale Numeric scalar. Scaling factor applied to the standard deviation
#'   of the normal distribution used for imputation. Default is \code{0.3}.
#' @param shift Numeric scalar. Left-shift factor applied to the median of the
#'   observed values (in units of standard deviation). Default is \code{1.8}.
#'
#' @return
#' A \code{proteoDA::DAList} object with missing values imputed in the
#' \code{data} slot.
#'
#' @importFrom stats cor na.omit quantile rnorm sd setNames
#'
#' @examples
#' \dontrun{
#' ## Impute missing values using default parameters
#' da_imputed <- mod_glob_imput_transform_fct_imput(
#'   DAlist_obj = da_list
#' )
#'
#' ## Custom imputation parameters
#' da_imputed <- mod_glob_imput_transform_fct_imput(
#'   DAlist_obj = da_list,
#'   scale = 0.25,
#'   shift = 2
#' )
#' }
#'
#' @export
mod_glob_imput_transform_fct_imput <- function(DAlist_obj, scale = 0.3, shift = 1.8) {
  scale <- as.numeric(scale)
  shift <- as.numeric(shift)
  # Show error if inputs are not the required class DALists
  assertthat::assert_that(inherits(DAlist_obj, "DAList"),
                          is.numeric(scale),
                          length(scale) == 1,
                          is.numeric(shift),
                          length(shift) == 1)

  DAList_assay <- DAlist_obj$data

  # Show error if there are no missing values
  if(!any(is.na(DAList_assay))) {
    stop("No missing values in '", deparse(substitute(DAlist_obj)), "'",
         call. = FALSE)
  }

  # Get descriptive parameters of the current sample distributions
  stat <- DAList_assay %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    gather(samples, value, -rowname) %>%
    filter(!is.na(value))  %>%
    group_by(samples) %>%
    summarise(mean = mean(value),
              median = median(value),
              sd = sd(value),
              n = n(),
              infin = nrow(DAList_assay) - n)
  # Impute missing values by random draws from a distribution
  # which is left-shifted by parameter 'shift' * sd and scaled by parameter 'scale' * sd.
  for (a in seq_len(nrow(stat))) {
    DAList_assay[is.na(DAList_assay[, stat$samples[a]]), stat$samples[a]] <-
      rnorm(stat$infin[a],
            mean = stat$median[a] - shift * stat$sd[a],
            sd = stat$sd[a] * scale)
  }
  DAlist_obj$data <- DAList_assay
  return(DAlist_obj)
}


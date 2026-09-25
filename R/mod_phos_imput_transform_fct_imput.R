#' Impute missing values in a DAList object
#'
#' This function imputes missing values in a phosphosite DAList object using a
#' left-shifted Gaussian distribution. The missing values are replaced by random
#' draws from a normal distribution with mean shifted below the observed median
#' and standard deviation scaled according to user-defined parameters.
#'
#' @param DAlist_obj A \code{DAList} object from \pkg{proteoDA} containing the
#'   phosphosite intensity data in the \code{$data} slot.
#' @param scale Numeric. Factor to scale the standard deviation of the left-shifted
#'   normal distribution. Default is 0.3.
#' @param shift Numeric. Factor to shift the mean of the normal distribution below
#'   the median of the observed values. Default is 1.8.
#'
#' @details
#' The imputation is performed per sample (column) as follows:
#' \enumerate{
#'   \item For each sample, compute the median and standard deviation of observed values.
#'   \item Count the number of missing values in that sample.
#'   \item Impute missing values by random draws from a normal distribution with:
#'     \itemize{
#'       \item mean = median - shift * sd
#'       \item sd   = scale * sd
#'     }
#' }
#'
#' This approach preserves the overall distribution of each sample while imputing
#' values likely below the detection limit.
#'
#' @return A \code{DAList} object with missing values in the \code{$data} slot replaced
#'   by imputed values.
#'
#' @importFrom stats rnorm sd
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import assertthat
#' @export
#'
#' @examples
#' \dontrun{
#' # Impute missing values in a DAList object
#' DAList_imputed <- mod_phos_imput_transform_fct_imput(
#'   DAlist_obj = my_DAList,
#'   scale = 0.3,
#'   shift = 1.8
#' )
#' }
mod_phos_imput_transform_fct_imput <- function(DAlist_obj, scale = 0.3, shift = 1.8) {
  if(is.integer(scale)) scale <- is.numeric(scale)
  if(is.integer(shift)) shift <- is.numeric(shift)
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


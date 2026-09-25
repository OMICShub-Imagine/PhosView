#' Calculate per-row variances of a numeric array
#'
#' Have done testing, gives same results as the genefilter::rowVars() function.
#'
#' @param x The array for which to calculate per-row variances
#' @param ... Additional arguments to be passed to internal functions.
#'   Meant for na.rm
#'
#' @return A numeric vector of appropriate length, named if input was named, with
#'   per-row variances.
#'
#' @keywords internal
#'
#' @examples
#' # No examples yet.

rowVars <- function(x, ...) {
  # Get variance of each row by vectorized subtraction of rowmeans
  # This is much faster than, e.g., applying the built-in var() function over the rows of a matrix
  rowSums((x - rowMeans(x, ...))^2, ...)/(rowSums(!is.na(x)) - 1)
}


#' Calculate per-row standard deviations of a numeric array
#'
#' Have done testing, gives same results as the genefilter::rowSds() function.
#'
#' @param x The array for which to calculate per-row Sds
#' @param ... Additional arguments to be passed to internal functions.
#'   Meant for na.rm.
#'
#' @return A numeric vector of appropriate length, named if input was named, with
#'   per-row standard deviations
#'
#' @keywords internal
#'
#' @examples
#' # No examples yet.
rowSds <- function(x, ...) {
  sqrt(rowVars(x, ...))
}

#' Calculate per-row medians of a numeric array
#'
#' Have done testing, gives same results as the matrixStats::rowMedians() function
#' it replaces, though it is much slower (the matrixStats version uses C code).
#'
#' @param x The array for which to calculate per-row medians
#' @param ... Additional arguments to be passed to internal functions.
#'   Meant for na.rm.
#'
#' @return A numeric vector of appropriate length, named if input was named, with
#'   per-row medians
#'
#' @keywords internal
#'
#' @examples
#' # No examples yet.
rowMedians <- function(x, ...) { # Much slower than the matrixStats version, which uses C
  # but probably worth losing the dependency.
  apply(x, MARGIN = 1, FUN = stats::median, ...)
}

#' Calculate per-row MADs of a numeric array
#'
#' Have done testing, gives same results as the matrixStats::rowMads() function
#' it replaces, though it is much slower (the matrixStats version uses C code).
#'
#' @param x The array for which to calculate per-row MADs
#' @param ... Additional arguments to be passed to internal functions.
#'   Meant for na.rm.
#'
#' @keywords internal
#'
#' @return A numeric vector of appropriate length, named if input was named, with
#'   per-row MADs
#'
#' @examples
#' # No examples yet.
rowMads <- function(x, ...) {
  apply(x, MARGIN = 1, FUN = stats::mad, ...)
}


#' Metrics to evaluate normalization methods
#'
#' A set of functions that take in normalized sample data and a grouping factor
#' and calculate some metric of of variability, error, etc., that we can use
#' to evaluate normalization methods. The functions/metrics: \itemize{
#'   \item PCV- Calculated the pooled coefficient of variation (CV, standard
#'     deviation divided by the mean) for all proteins within a group. That is,
#'     the average (mean) of the per-protein CVs.
#'   \item PMAD- Calculate the pooled median absolute deviation. The MAD is the
#'     median of the absolute deviations of each sample from the group median.
#'     The pooled MAD is the group mean of the per-protein PMADs.
#'   \item PEV- Calculates the pooled estimate of variance for each group. Uses
#'     the weighted average method to account for unequal sample sizes (see
#'     \url{https://en.wikipedia.org/wiki/Pooled_variance}).
#'   \item COR- Calculates all within-group pairwise correlations in intensity
#'     between samples.
#'   \item log2ratio- Calculates all possible between-group ratios of average
#'     log2-normalized intensity for each protein.
#' }
#'
#' @param data A numeric data frame, where each row is a protein and columns
#'   are densities. For most uses, these are probably normalized data.
#' @param groups A character or factor vector, listing the group(s) the samples
#'   belong to.
#'
#' @return For PCV, PMAD, and PEV: a named vector, with length equal to the
#'   number of groups, giving the metric for each group. For COR, a named vector
#'   giving all within-group pairwise correlations between samples. For log2ratio
#'   with keep_protein_ID = FALSE (the default), a non-named vector giving all
#'   possible between-group log2 ratios for each protein. For log2ratio with
#'   keep_protein_ID = TRUE, a data frame in which rows are proteins, columns are
#'   are the log2ratio for a pairwise comparisons, and rownames give protein ID.
#'
#' @name norm_metrics
#'
#' @rdname norm_metrics
#' @keywords internal
#'
# Pooled coefficient of variation (average of per-protein CVs)
PCV <- function(data, groups) {
  PCV <- NULL
  for (group in unique(groups)) {
    tempData <- as.matrix(data[, groups %in% group])
    CVs <- rowSds(tempData, na.rm = FALSE)/rowMeans(tempData, na.rm = FALSE)
    PCV[group] <- mean(CVs, na.rm = T)
  }
  return(PCV)
}

#' @rdname norm_metrics
#' @keywords internal
#'
# Pooled median absolute deviation (median of per-protein MADs)
PMAD <- function(data, groups) {
  PMAD <- NULL
  for (group in unique(groups)) {
    tempData <- as.matrix(data[, groups %in% group])
    MAD <- rowMads(tempData, na.rm = F)
    PMAD[group] <- mean(MAD, na.rm = T)
  }
  return(PMAD)
}


#' @rdname norm_metrics
#' @keywords internal
#'
# Pooled estimate of variance
PEV <- function(data, groups) {
  PEV <- NULL
  for (group in unique(groups)) {
    tempData <- as.matrix(data[, groups %in% group])

    rowNonNACnt <- rowSums(!is.na(tempData)) - 1
    EV <- rowNonNACnt * rowVars(tempData, na.rm = FALSE)
    PEV[group] <- sum(EV, na.rm = TRUE)/sum(rowNonNACnt, na.rm = TRUE)
  }
  return(PEV)

}

#' @rdname norm_metrics
#' @keywords internal
#'
# Within-group pairwise correlations.
COR <- function(data, groups) {
  COR <- NULL
  for (group in unique(groups)) { # for each group
    corGroup <- NULL
    tempData <- as.matrix(data[,groups %in% group])
    if (ncol(tempData) == 1) { # If only one sample, cor with itself is 1
      corGroup <- 1
    }
    if (ncol(tempData) > 1) { # If multiple sample, get all pw cors
      corVals <- stats::cor(tempData,
                            use = "pairwise.complete.obs",
                            method = "pearson")
      corGroup <- corVals[lower.tri(corVals)]
    }
    COR[[group]] <- corGroup
  }
  return(unlist(COR))
}

#' All pairwise differences of a vector
#'
#' Taken from \url{https://stackoverflow.com/questions/48445003/compute-all-pairwise-differences-of-elements-in-a-vector},
#' uses a nice R matrix algebra trick, taking the outer product of the vectors
#' (in this case, the difference instead of the product). Then, just return the
#' lower triangle, so we don't double-count pw comparisons
#'
#' @param vector A numeric vector of items for which you want all pairwise
#'   differences
#'
#' @return A vector of all pairwise differences between the elements in the vector
#'
#' @keywords internal
#'
#' @examples
#' # No examples yet
#'
all_pw_diffs <- function(vector) {
  all_diffs <- outer(vector, vector, "-")
  all_diffs[lower.tri(all_diffs)]
}



#' @rdname norm_metrics
#' @param keep_protein_ID Should protein ID information be retained? Default is FALSE.
#' @keywords internal
#'
# log2ratio
log2ratio <- function(data, groups, keep_protein_ID = F) {

  # calculate mean intensity for each group
  # Original implementation by charity was a nested for loop that ended up calculating group
  # means multiple times
  # I tried a one-liner with aggregate, but it was slow as hell
  # Switch back to a loop, but now a non-nested one that uses the speedy rowMeans
  mean_prot_by_group <- as.data.frame(matrix(nrow = nrow(data),
                                             ncol = length(unique(groups)),
                                             dimnames = list(rownames(data),
                                                             sort(unique(groups)))))
  # Get group means in a loop
  for (group in sort(unique(groups))) {
    one_group_data <- data[, groups == group, drop = F]
    mean_prot_by_group[, group] <- rowMeans(one_group_data, na.rm = T)
  }



  # Then, calculate all pairwise differences across cols/groups
  # Assuming these are all log2 normalized, the difference is the ratio
  # uses a custom function, all_pw_diffs (in utils), to apply this.
  # Need simplify = F and some ugly type conversion and transposition
  # To get this in the format we want (rows = proteins, cols = pw comparisons)
  # If simplify = T, cases with only 1 pw comparison are simplified to vectors.
  pw_diffs_per_gene <- as.data.frame(t(as.data.frame(apply(X = mean_prot_by_group, MARGIN = 1, FUN = all_pw_diffs, simplify = F))))


  # The result is a data frame, where each row is a gene and each col is the difference in one
  # pw comparison.
  # For the log2ratio plots, want to collapse protein info and just return
  # a vector like the other plots.
  # For the MD/MA plots, need to maintain protein info in a data frame
  if (keep_protein_ID) {
    output <- pw_diffs_per_gene
  } else {
    output <- unlist(pw_diffs_per_gene, use.names = F)
  }

  output
}


#' Evaluate a normalization metric across a list of normalized data
#'
#' Applies a metric (see \code{\link{norm_metrics}}) across a list of normalized
#' data matrices (most likely produced by \code{\link{apply_all_normalizations}})
#' and outputs the results in a data frame ready to be used with internal plotting
#' functions (see \code{\link{pn_plots_generic}}).
#'
#' @param normList A named list of normalized data matrices.
#' @param grouping A character or factor vector, listing the group(s) the samples
#'   belong to.
#' @param metric The normalization metric to calculate. Can be "PCV",
#'   "PMAD", "PEV", or "COR". See \code{\link{norm_metrics}}.
#' @return A data frame containing the selected metric,
#'   ready to be used for plotting.
#'
#' @keywords internal
#'
eval_pn_metric_for_plot <- function(normList,
                                    grouping,
                                    metric = c("PCV", "PMAD", "PEV", "COR", "log2ratio")) {

  # check args
  metric <- rlang::arg_match(metric)

  plotData <- NULL
  for (norm_method in names(normList)) {
    metric_results <- do.call(what = metric, args= list(data = as.data.frame(normList[norm_method]), groups = grouping))
    if (metric == "log2ratio") { #log2ratio does not have the per-group info the other metrics do, only per-normalization method
      x <- data.frame(method = norm_method,
                      value = metric_results, row.names = NULL)
    } else {
      x <- data.frame(method = norm_method,
                      group = names(metric_results),
                      value = metric_results, row.names = NULL)
    }
    plotData <- rbind(plotData, x)
  }

  plotData$method <- factor(plotData$method, levels = names(normList))

  return(plotData)
}


#' Generic plotting functions for normalization report
#'
#' An internal set of functions that take in processed data ready for plotting and
#' create generic plots (e.g., means, violin plots, line plots). These generic plots
#' are then passed to metric-specific plotting functions (see \code{\link{pn_plots}})
#' to create the final plot objects.
#'
#' @param plotData A data frame of data to be plotted, created with the
#'   \code{\link{eval_pn_metric_for_plot}} function.
#' @return A ggplot object of the plot.
#'
#' @name pn_plots_generic
#'
#' @importFrom ggplot2 ggplot aes
#' @importFrom ggplot2 geom_pointrange geom_violin geom_vline geom_line geom_point geom_hline geom_smooth
#' @importFrom ggplot2 position_nudge position_jitterdodge
#' @importFrom ggplot2 scale_fill_manual scale_color_manual
#' @importFrom ggplot2 theme_bw theme element_text element_blank element_rect
#' @importFrom ggplot2 xlab ylab ggtitle coord_cartesian unit facet_wrap
#'
#' @keywords internal
#' @rdname pn_plots_generic
#'
pn_mean_plot <- function(plotData) {

  blueberry <- "#1F5EDC";  cherry <- "#EE0010"
  apple     <- "#32CD32";  barbie <- "#FF1493"
  fanta     <- "#FF7F00";  grape  <- "#A342FC"
  ocean     <- "#00C8FF";  mtndew <- "#ADFF2F"
  gold      <- "#FFE100";  orchid <- "#E36EF6"
  aceblue   <- "#009ACE";  poop   <- "#996633"

  binfcolors <-c(blueberry, cherry, apple, barbie, fanta,
                 grape, ocean, mtndew, gold, orchid, aceblue, poop)
  names(binfcolors)<-c("blueberry", "cherry", "apple", "barbie", "fanta",
                       "grape", "ocean", "mtndew", "gold", "orchid", "aceblue", "poop")


  means <- stats::aggregate(value ~ method, data = plotData, mean)
  sds <- stats::aggregate(value ~ method, data = plotData, stats::sd)
  ns <- stats::aggregate(value ~ method, data = plotData, FUN = length)

  summary <- data.frame(method = means$method,
                        mean = means$value,
                        se = sds$value/sqrt(ns$value))

  # make the plot
  result <- summary  |>
    ggplot(aes(col = .data$method)) +
    geom_pointrange(aes(x = .data$method, y = .data$mean,
                        ymin = .data$mean - .data$se,
                        ymax = .data$mean + .data$se),
                    pch = 18,
                    size = 1.15,
                    position = position_nudge(x = -0.2)) +
    geom_point(aes(x = .data$method, .data$value),
               position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.3),
               data = plotData) +
    scale_color_manual(values = unname(binfcolors)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          panel.grid = element_blank())

  return(result)
}

#' @rdname pn_plots_generic
#'
pn_violin_plot <- function(plotData) {

  blueberry <- "#1F5EDC";  cherry <- "#EE0010"
  apple     <- "#32CD32";  barbie <- "#FF1493"
  fanta     <- "#FF7F00";  grape  <- "#A342FC"
  ocean     <- "#00C8FF";  mtndew <- "#ADFF2F"
  gold      <- "#FFE100";  orchid <- "#E36EF6"
  aceblue   <- "#009ACE";  poop   <- "#996633"

  binfcolors <-c(blueberry, cherry, apple, barbie, fanta,
                 grape, ocean, mtndew, gold, orchid, aceblue, poop)
  names(binfcolors)<-c("blueberry", "cherry", "apple", "barbie", "fanta",
                       "grape", "ocean", "mtndew", "gold", "orchid", "aceblue", "poop")

  # make the plot
  result <- plotData  |>
    ggplot(aes(fill = .data$method)) +
    geom_violin(aes(x = .data$method, y = .data$value),
                draw_quantiles = c(0.5),
                col = "black") +
    scale_fill_manual(values = unname(binfcolors)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          panel.grid = element_blank())

  return(result)
}

#' @rdname pn_plots_generic
#'
pn_density_plot <- function(plotData) {
  blueberry <- "#1F5EDC";  cherry <- "#EE0010"
  apple     <- "#32CD32";  barbie <- "#FF1493"
  fanta     <- "#FF7F00";  grape  <- "#A342FC"
  ocean     <- "#00C8FF";  mtndew <- "#ADFF2F"
  gold      <- "#FFE100";  orchid <- "#E36EF6"
  aceblue   <- "#009ACE";  poop   <- "#996633"

  binfcolors <-c(blueberry, cherry, apple, barbie, fanta,
                 grape, ocean, mtndew, gold, orchid, aceblue, poop)
  names(binfcolors)<-c("blueberry", "cherry", "apple", "barbie", "fanta",
                       "grape", "ocean", "mtndew", "gold", "orchid", "aceblue", "poop")
  # make the plot
  result <- plotData  |>
    ggplot(aes(color = .data$method, group = .data$method)) +
    geom_vline(aes(xintercept = 0), col = "grey80") +
    geom_line(aes(x = .data$value),
              stat = "density",
              na.rm = T,
              trim = T) + # need trim = T for the density calculation to happen per-group
    scale_color_manual(values = unname(binfcolors)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank())

  return(result)
}


#' Functions for plotting normalization metrics
#'
#' A set of functions that take in normalized sample data and a grouping factor
#' and calculate some metric of of variability, error, etc., that we can use
#' to evaluate normalization methods. See \code{\link{norm_metrics}} for explanations
#' of each metric.
#'
#' @inheritParams eval_pn_metric_for_plot
#' @param zoom Should the plot cover the full range of log2ratios, or zoom
#'   in around 0? Default is FALSE.
#' @param legend Should the plot include the legend? Default is TRUE.
#' @return A ggplot object of the plot.
#'
#' @name pn_plots
#'

#' @rdname pn_plots
#' @keywords internal
#'
pn_plot_PCV <- function(normList, grouping) {
  eval_pn_metric_for_plot(normList,
                          grouping,
                          metric = "PCV") |>
    pn_mean_plot() +
    ylab("Pooled Coefficient of Variation") +
    ggtitle("PCV")
}

#' @rdname pn_plots
#' @keywords internal
#'
pn_plot_PMAD <- function(normList, grouping) {
  eval_pn_metric_for_plot(normList,
                          grouping,
                          metric = "PMAD") |>
    pn_mean_plot() +
    ylab("Median Absolute Deviation") +
    ggtitle("PMAD")
}

#' @rdname pn_plots
#' @keywords internal
#'
pn_plot_PEV <- function(normList, grouping) {
  eval_pn_metric_for_plot(normList,
                          grouping,
                          metric = "PEV") |>
    pn_mean_plot() +
    ylab("Pooled Estimate of Variance") +
    ggtitle("PEV")
}

#' @rdname pn_plots
#' @keywords internal
#'
pn_plot_COR <- function(normList, grouping) {
  eval_pn_metric_for_plot(normList,
                          grouping,
                          metric = "COR") |>
    pn_violin_plot() +
    ylab("Intragroup Correlation") +
    ggtitle("COR")
}

#' @rdname pn_plots
#' @keywords internal
#'
pn_plot_log2ratio <- function(normList, grouping, zoom = F, legend = T) {

  # Get plot data
  plotData <- eval_pn_metric_for_plot(normList,
                                      grouping,
                                      metric = "log2ratio")

  # Find max density
  # and min/max of ratios
  # for zooming later
  maxY <- max(stats::aggregate(value ~ method, data = plotData, FUN = function(x) max(stats::density(x, na.rm = T)$y))$value)
  min_xlim <- 0.5 * min(stats::aggregate(value ~ method, data = plotData, FUN = function(x) min(stats::density(x, na.rm = T)$x))$value)
  max_xlim <- 0.5 * max(stats::aggregate(value ~ method, data = plotData, FUN = function(x) max(stats::density(x, na.rm = T)$x))$value)

  # Build base plot
  base <- plotData |>
    pn_density_plot() +
    xlab("") +
    ylab("Density") +
    ggtitle("Log2 ratio")

  # Do zooming
  if (zoom) {
    base <- base +
      coord_cartesian(xlim = c(-0.3, 0.3),
                      ylim = c(maxY - (0.5*maxY), maxY + (0.2*maxY)))
  } else {
    base <- base +
      coord_cartesian(xlim = c(min_xlim, max_xlim))
  }

  # Output baseplot with or without legend
  if (legend) {
    result <- base +
      theme(legend.justification = c(1,1), # sets upper-right corner as locating point
            legend.position = c(0.99, 0.99), # puts locating point in upper-right corner
            legend.title = element_blank(),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.35, "cm"),
            legend.background = element_rect(fill = NULL))
  } else {
    result <- base +
      theme(legend.position = "none")
  }

  result
}


#' @rdname pn_plots
#' @keywords internal
#'
pn_plot_MD <- function(normList, grouping, use_ggrastr = F) {

  # Assemble data for plotting
  log2ratios <- NULL
  for (norm_method in names(normList)) {
    metric_results <- do.call(what = "log2ratio", args= list(data = as.data.frame(normList[norm_method]),
                                                             groups = grouping, keep_protein_ID = T))

    long <- utils::stack(metric_results)
    long$protein  <- rep(rownames(metric_results), times = ncol(metric_results))
    long$method <- norm_method

    log2ratios <- rbind(log2ratios, long)
  }

  # Merge together the log2ratios for each pairwise comparison
  # and their average intensities
  plotData <- merge(log2ratios,
                    data.frame(mean_intensity = rowMeans(normList$log2, na.rm = T),
                               protein = rownames(normList$log2)))
  plotData$method <- factor(plotData$method, levels = names(normList))

  # make plot
  # possibly using ggrastr
  if (use_ggrastr) {
    if (!requireNamespace("ggrastr", quietly = TRUE)) {
      cli::cli_abort(c("Package \"ggrastr\" must be installed in order to use it.")) #nocov
    }

    base_plot <- plotData |>
      ggplot(aes(x = .data$mean_intensity, y = .data$values)) +
      ggrastr::geom_point_rast(alpha = 0.3, na.rm = T)
  } else {
    base_plot <- plotData |>
      ggplot(aes(x = .data$mean_intensity, y = .data$values)) +
      geom_point(alpha = 0.3, na.rm = T)

  }

  base_plot +
    geom_hline(yintercept = 0, color = "dodgerblue3") +
    geom_smooth(method = "lm", formula = "y ~ x", na.rm = T, color = "darkorange") +
    facet_wrap("method", ncol = 4) +
    ylab("logFC") +
    xlab("mean intensity (log2)") +
    theme_bw() +
    theme(panel.grid = element_blank())
}


#' DAList validator
#'
#' Internal function for validating a DAList object
#'
#' @param x An object to be tested if it is a valid DAList.
#'
#' @return If x was a valid DAList, returns x
#'
#' @keywords internal
#'
validate_DAList <- function(x) {

  ## Check overall structure

  slots <- c("data", "annotation", "metadata", "design",
             "eBayes_fit", "results", "tags")

  if (!identical(names(x), slots)) {

    missing_slots <- slots[!(slots %in% names(x))]
    if (length(missing_slots) > 0) {
      cli::cli_abort(
        "The DAList is missing the following slot{?s}: {missing_slots}"
      )
    }

    extra_slots <- names(x)[!(names(x) %in% slots)]
    if (length(extra_slots) > 0) {
      cli::cli_abort(
        "The DAList contains the following extra slot{?s}: {extra_slots}"
      )
    }

    cli::cli_alert("Slots in DAList out of order. Reordering.")
    x <- x[slots]
  }

  ## Check each element

  if (any(c("tbl_df", "tbl") %in% class(x$data))) {
    cli::cli_abort(
      "Data slot is a tibble, must be a matrix or data frame"
    )
  }

  if (any(c("tbl_df", "tbl") %in% class(x$annotation))) {
    cli::cli_abort(
      "Annotation slot is a tibble, must be a matrix or data frame"
    )
  }

  if (any(c("tbl_df", "tbl") %in% class(x$metadata))) {
    cli::cli_abort(
      "Metadata slot is a tibble, must be a matrix or data frame"
    )
  }

  if (!any(c(is.data.frame(x$data), is.matrix(x$data)))) {
    cli::cli_abort(
      "The {.arg data} slot of a DAList must be a data frame or matrix"
    )
  }

  if (!all(apply(x$data, 2, is.numeric))) {
    cli::cli_abort(
      "The {.arg data} slot of a DAList must contain only numeric data"
    )
  }

  ## Annotation checks

  if (!("uniprot_id" %in% colnames(x$annotation))) {
    cli::cli_abort(
      "The annotation data must contain a column named \"uniprot_id\""
    )
  }

  if (any(duplicated(x$annotation$uniprot_id))) {
    cli::cli_abort(c(
      "Entries in the {.val uniprot_id} column of the annotation are not unique.",
      "i" = "Find duplicate items with {.fun base::duplicated} or {.fun base::anyDuplicated} and make them unique"
    ))
  }

  if (nrow(x$data) != nrow(x$annotation)) {
    cli::cli_abort(
      "The {.arg data} slot and {.arg annotation} slots of a DAList must have the same number of rows"
    )
  }

  if (!all(rownames(x$data) == rownames(x$annotation))) {
    cli::cli_abort(
      "Rownames for the {.arg data} and {.arg annotation} slots of the DAList must match"
    )
  }

  if (!all(rownames(x$data) == x$annotation$uniprot_id)) {
    cli::cli_abort(
      "Rownames for the {.arg data} and {.arg annotation} slots of the DAList must equal the {.val uniprot_id} column of the annotation."
    )
  }

  ## Metadata checks

  if (!is.null(x$metadata)) {
    if (nrow(x$metadata) != ncol(x$data)) {
      cli::cli_abort(
        "The number of samples in the metadata ({nrow(x$metadata)}) do not match the number of samples in the data ({ncol(x$data)})"
      )
    }

    if (any(colnames(x$data) != rownames(x$metadata))) {
      cli::cli_abort(
        "The row names of the metadata do not match the column names of the data"
      )
    }
  }

  ## Design checks

  if (!is.null(x$design)) {

    if (length(x$design) < 2) {
      cli::cli_abort(c(
        "The object in the design slot should be a list with length of at least 2:",
        "Current length: {length(x$design)}"
      ))
    }

    if (!all(c("design_matrix", "design_formula") %in% names(x$design))) {
      cli::cli_abort(
        "The list in the design slot must contain at least a design_formula and a design_matrix"
      )
    }

    if (!("assign" %in% names(attributes(x$design$design_matrix)))) {
      cli::cli_abort(
        "The design_matrix in the design slot should have an \"assign\" attribute"
      )
    }

    design_possible_names <- c(
      "design_formula", "design_matrix", "random_factor",
      "contrast_matrix", "contrast_vector"
    )

    if (any(!(names(x$design) %in% design_possible_names))) {
      problem_names <- names(x$design)[!(names(x$design) %in% design_possible_names)]
      cli::cli_abort(c(
        "The list in the design slot may contain the following elements:",
        "{.val {design_possible_names}}",
        "Problem {cli::qty(problem_names)} item{?s} in design list: {.val {problem_names}}"
      ))
    }

    if (nrow(x$design$design_matrix) != nrow(x$metadata)) {
      cli::cli_abort(
        "The number of samples in the design matrix ({nrow(x$design$design_matrix)}) does not match the number of samples in the metadata ({nrow(x$metadata)})"
      )
    }

    if (any(colnames(x$data) != rownames(x$design$design_matrix))) {
      cli::cli_abort(
        "The row names of the design matrix do not match the column names of the data"
      )
    }

    if (!is.null(x$design$random_factor)) {
      if (!(x$design$random_factor %in% colnames(x$metadata))) {
        cli::cli_abort(
          "The random factor term is not present in the metadata"
        )
      }
    }

    if (sum(c("contrast_vector", "contrast_matrix") %in% names(x$design)) > 0) {
      if (sum(c("contrast_vector", "contrast_matrix") %in% names(x$design)) != 2) {
        cli::cli_abort(
          "If contrast information is present in the design slot, both a contrast_vector and contrast_matrix must be provided"
        )
      }

      if (!all(rownames(x$design$contrast_matrix) ==
               colnames(x$design$design_matrix))) {
        cli::cli_abort(
          "The rows in the contrast matrix do not match the columns in the design matrix"
        )
      }
    }
  }

  ## eBayes checks

  if (!is.null(x$eBayes_fit)) {

    if (!("MArrayLM" %in% class(x$eBayes_fit))) {
      cli::cli_abort(
        "The eBayes slot should be an object of class MArrayLM"
      )
    }

    if (!is.null(x$design$random_factor)) {
      if (!("correlation" %in% names(x$eBayes_fit))) {
        cli::cli_abort(c(
          "Discrepancy between limma fit and statistical design:",
          "Design includes a random factor not present in the model fit"
        ))
      }
    }

    if ("correlation" %in% names(x$eBayes_fit)) {
      if (is.null(x$design$random_factor)) {
        cli::cli_abort(c(
          "Discrepancy between limma fit and statistical design:",
          "Model fit includes a random factor not present in the design"
        ))
      }
    }
  }

  ## Results checks

  if (!is.null(x$results)) {

    if (!is.null(x$design$contrast_matrix)) {
      expected_terms <- colnames(x$design$contrast_matrix)
      terms_from <- "contrast matrix"
    } else {
      expected_terms <- colnames(x$design$design_matrix)
      terms_from <- "design matrix"
    }

    if (x$tags$extract_intercept) {
      if (!all(names(x$results) == expected_terms)) {
        cli::cli_abort(c(
          "Discrepancy between statistical results and design:",
          "{cli::qty(length(expected_terms))} Term name{?s} in {terms_from} do not match {cli::qty(length(names(x$results)))} term{?s} names in the results"
        ))
      }
    } else {
      non_intercept_terms <- expected_terms[
        !stringr::str_detect(expected_terms, "ntercept")
      ]

      if (!all(names(x$results) == non_intercept_terms)) {
        cli::cli_abort(c(
          "Discrepancy between statistical results and design:",
          "{cli::qty(length(expected_terms))} Non-intercept term name{?s} in {terms_from} do not match {cli::qty(length(names(x$results)))} term{?s} names in the results"
        ))
      }
    }

    for (term in names(x$results)) {
      if (!all(rownames(x$results[[term]]) == rownames(x$data))) {
        cli::cli_abort(c(
          "Discrepancy between statistical results and data:",
          "Rownames for the {.arg {term}} results data frame do not match the rownames in the {.arg data}"
        ))
      }
    }
  }

  x
}


#' Create a normalization report
#' This is a modified version of the original normalization_report.R script from \pkg{proteoDA}
#' This modified version allows to generate the normalization report in a format compatible
#' with \pkg{shiny} rendering (no PDF creation but use of \pkg{plotly}. In addition, no need of
#' \code{output_dir}, \code{overwrite}, \code{filename}, \code{use_ggrastr} arguments
#'
#' @param DAList A DAList.
#' @param grouping_column The name of the column in the metadata which
#'   gives information on how to group samples for normalization. Must be
#'   supplied: some metrics can't be calculated for only one group.
#' @param suppress_zoom_legend Should the legend be removed from the zoomed
#'   log2ratio plot? Default is FALSE.
#'
#' @return If report is created successfully, invisibly returns the
#'   input DAList.
#'
#' @importFrom ggplot2 ggsave
#'
#' @export
#'

mod_phos_norm_report_fct_norm_report2 <- function(DAList,
                                                  grouping_column = NULL,
                                                  suppress_zoom_legend = FALSE) {



  #################################
  ## Check args and set defaults ##
  #################################

  input_DAList <- validate_DAList(DAList)

  if (!is.null(DAList$tags$normalized)) {
    if (DAList$tags$normalized) {
      cli::cli_abort("Data in DAList are already normalized. Cannot evaluate normalization metrics")
    }
  }

  # If provided, check that grouping column exists in the target data frame
  # And set it
  if (is.null(grouping_column)) { # If no groups provided, abort
    cli::cli_abort("{.arg grouping_column} cannot be empty")
  } else { # check that grouping column exists in the target data frame
    if (length(grouping_column) != 1) {
      cli::cli_abort(c("Length of {.arg grouping_column} does not equal 1",
                       "i" = "Only specify one column name for {.arg grouping_column}"))

    }
    if (!grouping_column %in% colnames(DAList$metadata)) {
      cli::cli_abort(c("Column {.arg {grouping_column}} not found in metadata of {.arg DAList}",
                       "i" = "Check the column names with {.code colnames(DAList$metadata)}."))
    }
    # And set it
    groups <- as.character(DAList$metadata[,grouping_column])
    # And give error if there's only one group
    if (length(unique(groups)) < 2) {
      cli::cli_abort(c("Column {.arg {grouping_column}} does not contain at least two different groups",
                       "!" = "Cannot calculate all normalization metrics without at least two groups"))
    }
  }


  ###########################
  ## Do all normalizations ##
  ###########################
  cli::cli_inform("Starting normalizations")
  normList <- apply_all_normalizations(DAList$data)
  cli::cli_inform("Normalizations finished")


  ###########################
  ## Make the plot objects ##
  ###########################

  # First page: all the combined plotting metrics
  list(
    a = pn_plot_PCV(normList, groups),
    b = pn_plot_PMAD(normList, groups),
    c = pn_plot_PEV(normList, groups),
    d = pn_plot_COR(normList, groups),
    e = pn_plot_log2ratio(normList, groups),
    f = pn_plot_log2ratio(normList, groups, zoom = T, legend = !suppress_zoom_legend),
    g = pn_plot_MD(normList, groups)
  )

}

#' @keywords internal
#'
log2Norm <- function(dat) {
  logInt <- log2(as.matrix(dat))
  logInt[is.infinite(as.matrix(logInt))] <- NA
  return(as.matrix(logInt))
}

#' @keywords internal
#'
medianNorm <- function(logDat) {
  # Find medians of each sample
  # Divide by median
  # Multiply by mean of medians
  sampleMed <- apply(logDat, 2, stats::median, na.rm = TRUE)
  meanMed <- mean(sampleMed, na.rm = TRUE)
  out <- t(t(logDat) / sampleMed)
  out <- out * meanMed
  return(as.matrix(out))
}

#' @keywords internal
#'
meanNorm <- function(logDat) {
  # Find means of each sample
  # Divide by mean
  # Multiply by mean of means
  sampleMean <- apply(logDat, 2, mean, na.rm = TRUE)
  meanMean <- mean(sampleMean, na.rm = TRUE)
  out <- t(t(logDat) / sampleMean)
  out <- out * meanMean
  return(as.matrix(out))
}

#' @keywords internal
#'
vsnNorm <- function(dat) {
  vsnNormed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
  colnames(vsnNormed) <- colnames(dat)
  row.names(vsnNormed) <- rownames(dat)
  return(as.matrix(vsnNormed))
}

#' @keywords internal
#'
quantileNorm <- function(logDat) {
  quantNormed <- preprocessCore::normalize.quantiles(as.matrix(logDat), copy = TRUE)
  colnames(quantNormed) <- colnames(logDat)
  row.names(quantNormed) <- rownames(logDat)
  return(as.matrix(quantNormed))
}

#' @keywords internal
#'
cycloessNorm <- function(logDat) {
  cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(logDat), method = "fast")
  colnames(cycLoessNormed) <- colnames(logDat)
  row.names(cycLoessNormed) <- rownames(logDat)
  return(as.matrix(cycLoessNormed))
}

#' @keywords internal
#'
rlrNorm <- function(logDat) {
  # First, an internal function to do get the coefficients from rlm regression
  # for one column
  get_rlm_coeffs_and_slopes <- function(unnormalized_col, predictors) {
    MASS::rlm(unnormalized_col ~ predictors,
              na.action = stats::na.exclude)$coefficients
  }

  # Median intensities across samples are the predictors for the rlm
  row_medians <- rowMedians(logDat, na.rm = TRUE)

  # Apply the one-col function across the input matrix cols
  coefficients <- t(apply(X = logDat,
                          MARGIN = 2,
                          FUN = get_rlm_coeffs_and_slopes,
                          predictors = row_medians))
  stopifnot(nrow(t(logDat)) == nrow(coefficients))
  # Do the normalization: subtract the intercept from each point, then divide
  # by slope
  normalized_mat <- t((t(logDat) - coefficients[,1])/coefficients[,2])
  normalized_mat
}

#' @keywords internal
#'
giNorm <- function(dat) {
  # Make sure input data is matrix
  dat_mat <- as.matrix(dat)
  # Follow same pattern as mean and median, but with
  # sums
  col_sums <- colSums(dat_mat, na.rm = TRUE)
  median_col_sum <- stats::median(col_sums)

  normMatrix <- t(t(dat_mat)/col_sums) * median_col_sum
  # extra step to return log2 values, as with previous implementation
  normLog2Matrix <- log2(normMatrix)
  normLog2Matrix[is.infinite(normLog2Matrix)] <- NA
  normLog2Matrix
}


#' Apply all normalization methods to a set of raw data
#'
#' Takes in raw intensity data and applies 8 different normalization methods,
#' returning a list of the data normalized with each method.
#'
#' @param data A data frame or matrix of raw data to be normalized. Rows are proteins and
#'   columns are raw intensity data.
#'
#' @return A list length 8, where each item in the list is a
#'       named matrix. Names give the normalization method, and the matrix
#'       contains the normalized data.
#'
#' @keywords internal
#'
apply_all_normalizations <- function(data) {

  normList <- NULL

  ## apply normalization methods using functions listed above.
  ## NOTE: most of the normalizations use log2(intensity) as input except
  ## VSN which normalizes using raw intensities with no log2 transformation.
  normList[["log2"]]     <- log2Norm(dat = data)
  normList[["median"]]   <- medianNorm(logDat = normList[["log2"]])
  normList[["mean"]]     <- meanNorm(logDat = normList[["log2"]])
  normList[["vsn"]]      <- vsnNorm(dat = data)
  normList[["quantile"]] <- quantileNorm(logDat = normList[["log2"]])
  normList[["cycloess"]] <- cycloessNorm(logDat = normList[["log2"]])
  normList[["rlr"]]      <- rlrNorm(logDat = normList[["log2"]])
  normList[["gi"]]       <- giNorm(dat = data)

  normList
}

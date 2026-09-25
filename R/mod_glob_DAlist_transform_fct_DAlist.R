#' Transform quantitative data and metadata into a proteoDA DAList object
#'
#' @description
#' Prepares quantitative protein group data and experimental metadata for
#' differential abundance analysis by constructing a \code{proteoDA::DAList}
#' object. The function extracts intensity values, cleans invalid entries,
#' validates sample naming consistency, and builds the required intensity,
#' annotation, and metadata components.
#'
#' @details
#' The transformation workflow includes:
#' \itemize{
#'   \item Determining the number of samples from the metadata.
#'   \item Extracting quantitative intensity columns from the input data.
#'   \item Replacing \code{NaN}, values < 1, and \code{"Filtered"} entries with
#'   \code{NA}.
#'   \item Creating an annotation table from protein group-level metadata.
#'   \item Duplicating \code{PG.ProteinGroups} into a mandatory
#'   \code{uniprot_id} column required by \code{proteoDA}.
#'   \item Verifying that sample names and ordering match between quantitative
#'   data and metadata.
#'   \item Renaming intensity columns using simplified sample identifiers.
#' }
#'
#' Summary statistics of the quantitative data are printed before and after
#' transformation to facilitate quality control.
#'
#' @param input A \code{data.frame} containing quantitative protein group data,
#'   typically loaded using
#'   \code{\link{mod_glob_DAlist_browser_fct_load}} with
#'   \code{data_type = "quant"}.
#' @param sample_metadata A \code{data.frame} containing experimental metadata.
#'   It must include at least the columns \code{Sample_ID} and
#'   \code{Simplify_ID}, with rows corresponding to samples.
#'
#' @return
#' A \code{proteoDA::DAList} object containing:
#' \itemize{
#'   \item \code{data}: cleaned quantitative intensity matrix
#'   \item \code{annotation}: protein group annotations
#'   \item \code{metadata}: experimental design information
#' }
#'
#' The function stops with an error if sample names or their ordering do not
#' match between quantitative data and metadata.
#'
#' @importFrom dplyr select mutate_all contains all_of
#' @importFrom tibble column_to_rownames
#' @importFrom skimr skim
#'
#' @examples
#' \dontrun{
#' ## Transform data into a DAList object
#' da_list <- mod_glob_DAlist_transform_fct_DAlist(
#'   input = quant_df,
#'   sample_metadata = meta_df
#' )
#' }
#'
#' @export


mod_glob_DAlist_transform_fct_DAlist <- function(input, sample_metadata) {

  ### Set a sample_number variable according to project based on exp_design.txt file number of rows
  sample_number <- dim(sample_metadata)[1]

  if ("PG.ProteinGroups" %in% colnames(input)) {

  #### Create intensity_data and annotation data
  intensity_data <- input %>%
    dplyr::select(c(11:(10+sample_number), PG.ProteinGroups)) %>%
    tibble::column_to_rownames("PG.ProteinGroups") %>%
    dplyr::mutate_all(~ifelse(is.nan(.), NA, .)) %>% ### replace NaNs by NA for further compatibility (normalisation step)
    dplyr::mutate_all(~ifelse(as.numeric(.) < 1, NA, .)) %>%  # Replace values <1 with NA
    dplyr::mutate_all(~ifelse(. == "Filtered", NA, .)) %>% ## replace "Filtered" cells with NA
    as.data.frame()

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

  annotation_data <- input %>%
    dplyr::select(all_of(required_cols)) %>%
    dplyr::mutate_all(~ifelse(is.nan(.), NA, .)) %>% ### replace NaNs by NA for further compatibility (normalisation step)
    as.data.frame()


  ### DAList needs to have a "uniprot_id" column name so duplicates "PG.ProteinGroups" into "uniprot_id"
  annotation_data$uniprot_id <- annotation_data$PG.ProteinGroups

  # Check if names and order match
  if (identical(colnames(intensity_data), sample_metadata$Sample_ID)) {
    print("The sample names and the order are the same between quantitative data and metadata.")
  } else {
    stop("The sample names or the order do not match between quantitative data and metadata, please check your Sample_ID")
  }

  ### Modified column name of intensity_data, in order to simplify them thanks to Simplify_ID column from exp_design
  colnames(intensity_data) <- sample_metadata$Simplify_ID

  # Check if names and order match
  if (identical(colnames(intensity_data), sample_metadata$Simplify_ID)) {
    print("The newly simplify IDs and the order are the same between quantitative data and metadata.")
  } else {
    stop("The newly simplify IDs or the order do not match between quantitative data and metadata, please check your Simplify_ID")
  }


  ### Put column name of quantitative value as rownames of sample_metadata object
  rownames(sample_metadata) <- sample_metadata$Simplify_ID

  } else if ("Protein.Group" %in% colnames(input)) {

    #### Create intensity_data and annotation data
    intensity_data <- input %>%
      dplyr::select(c(7:(6+sample_number), Protein.Group)) %>%
      tibble::column_to_rownames("Protein.Group") %>%
      dplyr::mutate_all(~ifelse(is.nan(.), NA, .)) %>% ### replace NaNs by NA for further compatibility (normalisation step)
      as.data.frame()

    required_cols <- c(
      "Protein.Group",
      "Protein.Names",
      "Genes",
      "First.Protein.Description",
      "N.Sequences",
      "N.Proteotypic.Sequences"
    )

    annotation_data <- input %>%
      dplyr::select(all_of(required_cols)) %>%
      dplyr::mutate_all(~ifelse(is.nan(.), NA, .)) %>% ### replace NaNs by NA for further compatibility (normalisation step)
      as.data.frame()


    ### DAList needs to have a "uniprot_id" column name so duplicates "PG.ProteinGroups" into "uniprot_id"
    annotation_data$uniprot_id <- annotation_data$Protein.Group

    # Check if names and order match
    if (identical(colnames(intensity_data), sample_metadata$Sample_ID)) {
      print("The sample names and the order are the same between quantitative data and metadata.")
    } else {
      stop("The sample names or the order do not match between quantitative data and metadata, please check your Sample_ID")
    }

    ### Modified column name of intensity_data, in order to simplify them thanks to Simplify_ID column from exp_design
    colnames(intensity_data) <- sample_metadata$Simplify_ID

    # Check if names and order match
    if (identical(colnames(intensity_data), sample_metadata$Simplify_ID)) {
      print("The newly simplify IDs and the order are the same between quantitative data and metadata.")
    } else {
      stop("The newly simplify IDs or the order do not match between quantitative data and metadata, please check your Simplify_ID")
    }


    ### Put column name of quantitative value as rownames of sample_metadata object
    rownames(sample_metadata) <- sample_metadata$Simplify_ID

  }

  raw <- proteoDA::DAList(data = intensity_data,
                annotation = annotation_data,
                metadata = sample_metadata,
                design = NULL,
                eBayes_fit = NULL,
                results = NULL,
                tags = NULL)

  return(raw)

}


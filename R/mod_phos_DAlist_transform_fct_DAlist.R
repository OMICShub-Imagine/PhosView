#' Build a proteoDA DAList object for phosphoproteomics data
#'
#' This helper function converts phosphosite quantitative data and its
#' experimental design into a \code{proteoDA::DAList} object. It performs
#' consistency checks between the phosphosite intensity matrix and the
#' experimental design, restructures quantitative and annotation data,
#' standardizes column names using simplified sample identifiers, and prepares
#' the data structure required for downstream differential analysis.
#'
#' @param quant_data A data.frame or tibble containing phosphosite-level
#'   quantitative data. Must include columns matching the pattern
#'   \code{"PTM.Quantity"} and a unique identifier column \code{PTM.CollapseKey}.
#'
#' @param meta_data A data.frame describing the experimental design. Must
#'   contain at least the columns \code{Sample_ID} and \code{Simplify_ID}, where
#'   \code{Sample_ID} corresponds to the original quantitative column order.
#'
#' @details
#' The function checks that the order of phosphosite quantitative columns
#' (matching \code{"PTM.Quantity"}) exactly matches the order of
#' \code{Sample_ID} in the experimental design. If this condition is not met,
#' execution is stopped to prevent sample misalignment.
#'
#' Quantitative columns are extracted and stored as the \code{data} slot of the
#' \code{DAList}, while remaining columns are used as annotation. The
#' phosphosite identifier \code{PTM.CollapseKey} is used as row names and
#' duplicated into a required \code{uniprot_id} annotation column.
#'
#' Sample column names are simplified using \code{Simplify_ID}, and PTM-related
#' annotation columns (site probability and PTM group) are renamed accordingly
#' to preserve sample-level consistency.
#'
#' @return A \code{proteoDA::DAList} object containing:
#' \itemize{
#'   \item \code{data}: phosphosite quantitative matrix (log-scale assumed)
#'   \item \code{annotation}: phosphosite metadata including PTM annotations
#'   \item \code{metadata}: experimental design information
#' }
#'
#' @seealso \code{\link[proteoDA]{DAList}}
#'
#' @export
mod_phos_DAlist_transform_fct_DAlist <- function(quant_data, meta_data) {

  ### Set a sample_number variable according to project based on exp_design.txt file number of rows
  sample_number <- dim(meta_data)[1]

  if ("Protein" %in% colnames(quant_data)) {

    #### Create intensity_data and annotation data
    intensity_data <- quant_data %>%
      dplyr::select(c(7:(6 + sample_number), Protein, Gene.Names, Residue, Site)) %>%
      dplyr::mutate(
        # Choose Gene.Names if available, otherwise fallback to Protein
        ID = dplyr::if_else(
          !is.na(Gene.Names) & Gene.Names != "",
          Gene.Names,
          Protein
        ),
        PTM.CollapseKey = paste0(ID, "_", Residue, Site)
      ) %>%
      dplyr::select(-c(Protein, Gene.Names, Residue, Site, ID)) %>%
      tibble::column_to_rownames("PTM.CollapseKey") %>%
      as.data.frame()

    intensity_data[intensity_data == 0] <- NA

    required_cols <- c(
      "Protein",
      "Protein.Names",
      "Gene.Names",
      "Residue",
      "Site",
      "Sequence"
    )

    annotation_data <- quant_data %>%
      dplyr::select(all_of(required_cols)) %>%
      dplyr::mutate(
        # Choose Gene.Names if available, otherwise fallback to Protein
        ID = dplyr::if_else(
          !is.na(Gene.Names) & Gene.Names != "",
          Gene.Names,
          Protein
        ),
        PTM.CollapseKey = paste0(ID, "_", Residue, Site)
      ) %>%
      dplyr::mutate_all(~ifelse(is.nan(.), NA, .)) %>% ### replace NaNs by NA for further compatibility (normalisation step)
      as.data.frame()


    ### DAList needs to have a "uniprot_id" column name so duplicates "PTM.CollapseKey" into "uniprot_id"
    annotation_data$uniprot_id <- annotation_data$PTM.CollapseKey

    # Check if names and order match
    if (identical(colnames(intensity_data), meta_data$Sample_ID)) {
      print("The sample names and the order are the same between quantitative data and metadata.")
    } else {
      stop("The sample names or the order do not match between quantitative data and metadata, please check your Sample_ID")
    }

    ### Modified column name of intensity_data, in order to simplify them thanks to Simplify_ID column from exp_design
    colnames(intensity_data) <- meta_data$Simplify_ID

    # Check if names and order match
    if (identical(colnames(intensity_data), meta_data$Simplify_ID)) {
      print("The newly simplify IDs and the order are the same between quantitative data and metadata.")
    } else {
      stop("The newly simplify IDs or the order do not match between quantitative data and metadata, please check your Simplify_ID")
    }


    ### Put column name of quantitative value as rownames of sample_metadata object
    rownames(meta_data) <- meta_data$Simplify_ID

  } else if ("PG.ProteinGroups" %in% colnames(quant_data)) {

  if(identical(colnames(quant_data)[grep("PTM.Quantity", colnames(quant_data))], meta_data$Sample_ID)) {
    print("The order of phosphosite matrix column names matches the order of Sample_ID in exp_design file.")
  } else {
    stop("The order of phosphosite matrix column names does not match the order of Sample_ID in exp_design file.")
  }

  ### Set a sample_number variable according to project based on exp_design.txt file number of rows

  #### Create intensity_data and annotation data and set a unique ID ("name") avec DEP::make_unique functionL
  quantity_cols <- grep("PTM\\.Quantity", names(quant_data), value = TRUE)

  intensity_data <- quant_data  %>%
    dplyr::select(all_of(c(quantity_cols, "PTM.CollapseKey"))) %>%
    column_to_rownames("PTM.CollapseKey")

  annotation_data <- quant_data  %>%
    dplyr::select(!all_of(quantity_cols)) %>%
    as.data.frame()

  ### DAList needs to have a "uniprot_id" column name so duplicates "PTM_collapse_key" into "uniprot_id"
  annotation_data$uniprot_id <- annotation_data$PTM.CollapseKey

  ### Modified column name of intensity_data, in order to simplify them thanks to Simplify_ID column from exp_design
  colnames(intensity_data) <- meta_data$Simplify_ID

  ### Put column name of quantitative value as rownames of meta_data object
  rownames(meta_data) <- meta_data$Simplify_ID

  ## Rename colnames of site prob and PTM.group as "Simplify_ID_PTM.SiteProbability" and "Simplify_ID_PTM.Group"
  siteprob_cols <- grep("PTM\\.SiteProbability", names(annotation_data))
  colnames(annotation_data)[siteprob_cols] <- paste0(meta_data$Simplify_ID, "_PTM_SiteProbability")

  ## Rename colnames of site prob and PTM.group as "Simplify_ID_PTM.SiteProbability" and "Simplify_ID_PTM.Group"
  PTMgroup_cols <- grep("PTM\\.Group", names(annotation_data))
  colnames(annotation_data)[PTMgroup_cols] <- paste0(meta_data$Simplify_ID, "_PTM_Group")

  }


  raw <- proteoDA::DAList(data = intensity_data,
                annotation = annotation_data,
                metadata = meta_data,
                design = NULL,
                eBayes_fit = NULL,
                results = NULL,
                tags = NULL)


  return(raw)

}


#' Prepare QC data for plotting
#'
#' @param quant_data should corresponds to a DAlist instance (DAlist_obj$data)
#' @param meta_data  should corresponds to a DAlist instance (DAlist_obj$metadata)
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom utils data head read.table str
#' @return A list of prepared data for QC plots

mod_phos_first_class_qc_utils_data_prep <- function(quant_data, meta_data) {

  ####################################################################################################################
  ###### TRANSFORM PHOSPHOSITE MATRIX (NUMERICAL VALUE RECOGNIZED AS CHARACTERS BECAUSE OF THE "FILTERING" TAG) ######
  ####################################################################################################################

  data <- quant_data
  metadata <- meta_data

  ### Set a sample_number variable according to project based on exp_design.txt file number of rows
  sample_number <- dim(metadata)[1]

  # Display structure before transformation
  cat("\n--- Structure BEFORE Transformation ---\n")
  str(data)

  # Reassign input to input_tf (tf = transformed)
  input_tf <- data

  # Replace "Filtered" with NA in all columns
  input_tf[input_tf == "Filtered"] <- NA

  cat("All Filtered tags were replace by NA.\n")

  # Function to replace "," with "." and convert to numeric
  convert_to_numeric <- function(x) {
    x <- gsub(",", ".", x)  # Replace comma with dot
    as.numeric(x)           # Convert to numeric
  }

  # Identify columns that contain "PTM.SiteProbability" or "PTM.Quantity"
  cols_to_convert <- grep("PTM\\.SiteProbability|PTM\\.Quantity", names(input_tf), value = TRUE)

  # Validate if columns exist
  if (length(cols_to_convert) == 2*sample_number) {
    cat("\nFound", length(cols_to_convert), "columns to transform:\n", paste(cols_to_convert, collapse = ",\n"), "\n")

    # Save original column classes for comparison
    original_classes <- sapply(input_tf[cols_to_convert], class)

    # Apply conversion to selected columns
    input_tf[cols_to_convert] <- lapply(input_tf[cols_to_convert], convert_to_numeric)

    # Check transformation results
    new_classes <- sapply(input_tf[cols_to_convert], class)

    # Compare before and after
    for (col in cols_to_convert) {
      if (original_classes[col] != "numeric" && new_classes[col] == "numeric") {
        cat("Column", col, "was successfully converted to numeric.\n")
      } else {
        cat("Column", col, "was NOT successfully converted.\n")
      }
    }
  } else {
    cat("\nNo matching columns found for conversion!\n")
  }

  ###########################################################################################################################
  ###### TRANSFORM PHOSPHOSITE MATRIX (ALL RAW INTENSITY < 1 TRANSFORM IN NA TO AVOID NEGATIVE VALUE AFTER LOG2 LATER) ######
  ###########################################################################################################################

  # Identify columns that contain "PTM.Quantity"
  quantity_cols <- grep("PTM\\.Quantity", names(input_tf), value = TRUE)

  # Replace numerocal value <1 with NA in PTM.Quantity columns only
  if (length(quantity_cols) == sample_number) {
    for (col in quantity_cols) {
      input_tf[[col]] <- ifelse(input_tf[[col]] < 1, NA, input_tf[[col]])
      cat("Column", col, "replaced quant value <1 by NA.\n")
    }
  } else {
    cat("No PTM.Quantity columns found to replace quant value <1 by NA.\n")
  }

  # Display structure after transformation
  cat("\n--- Structure AFTER Transformation ---\n")
  str(input_tf)
  summary(input_tf[quantity_cols])

  #################################################################
  ###### REMOVE PEPTIDE WITH OTHER MODIFICATION THAN PHOSPHO ######
  #################################################################

  ### Remove rows with PTM.Modification.Title different than "Phospho"

  # Create a new object input_tf_phospho as a copy of the original data
  input_tf_phospho <- input_tf

  # Find rows with modifications different from Phospho
  non_phospho_rows <- !grepl("Phospho", input_tf_phospho$PTM.ModificationTitle, ignore.case = TRUE)

  # Count the number of rows with non-Phospho modifications
  removed_rows <- sum(non_phospho_rows)

  # Remove rows with non-Phospho modifications
  input_tf_phospho <- input_tf_phospho[!non_phospho_rows, ]

  # Output message with the number of removed rows
  cat(removed_rows, "over",dim(input_tf)[1],"rows with modifications other than phosphorylation were removed.\n")

  ####################################################################################################
  ### Filtered rows where PTM.SiteProbability > 0.75 in X % of the replicate in at least one group ###
  ####################################################################################################

  # Create a new object First_class as a copy of input_tf
  first_class_tmp <- input_tf_phospho

  calculate_site_probability_stats <- function(first_class_tmp, metadata) {

    # Step 1: Identify PTM.SiteProbability columns
    site_prob_cols <- grep("PTM\\.SiteProbability", names(first_class_tmp), value = TRUE)

    # Step 2: Extract sample names from column names
    sample_names <- sub("\\.PTM\\.SiteProbability", "", site_prob_cols)

    # Step 4: Reshape first_class_tmp for processing
    site_prob_data <- first_class_tmp[, c("PTM.CollapseKey", site_prob_cols)]
    site_prob_data_long <- site_prob_data %>%
      pivot_longer(
        cols = -PTM.CollapseKey,  # Select all columns except PTM.CollapseKey
        names_to = "Sample_ID",
        values_to = "Site_Probability"
      )

    # Step 5: Merge with sample metadata
    site_prob_data_long$Sample_ID <- sub("\\.PTM\\.SiteProbability", "", site_prob_data_long$Sample_ID)
    metadata$Sample_ID <- sub("\\.PTM\\.Quantity", "", metadata$Sample_ID)
    merged_data <- left_join(site_prob_data_long, metadata, by = "Sample_ID")

    # Step 6: Calculate stats per group
    summary_stats <- merged_data %>%
      group_by(group, PTM.CollapseKey) %>%
      summarise(
        Non_Missing_Count = sum(!is.na(Site_Probability)),
        Count_LocProbability_Above_0.75 = sum(Site_Probability > 0.75, na.rm = TRUE)
      ) %>%
      ungroup()

    return(summary_stats)
  }

  result <- calculate_site_probability_stats(first_class_tmp, metadata)

  phosphosite_first_class <- function(result, metadata, replicate_threshold) {

    # Compute number of replicate per group
    rep_per_grp <- metadata %>%
      group_by(group) %>%
      summarise(rep_per_grps = n())

    # Map rep_per_group in data
    data <- left_join(result, rep_per_grp, by = "group")

    # Compute the minimum required count based on the threshold
    min_count <- round(data$rep_per_grps * (replicate_threshold / 100))

    # add min_count in data
    data <- cbind(result, min_count)

    # Compute number of replicate per group

    # Filter rows where at least one group meets the condition
    first_class_site <- data %>%
      filter(Non_Missing_Count >= min_count) %>%
      filter(Count_LocProbability_Above_0.75 >= min_count) %>%
      dplyr::select("PTM.CollapseKey") %>%
      unique() %>%
      unname() %>%
      unlist()

    return(first_class_site)
  }


  # Keep phosphosites with localization probability > 75% in at least 50% of replicates
  phospho_list_50 <- phosphosite_first_class(result, metadata, 50)
  phospho_list_70 <- phosphosite_first_class(result, metadata, 70)
  phospho_list_100 <- phosphosite_first_class(result, metadata, 100)

  ## Filter original dataframe based on IDs
  first_class_unfiltered <- first_class_tmp

  first_class_50 <- first_class_tmp %>%
    filter(PTM.CollapseKey %in% phospho_list_50)

  first_class_70 <- first_class_tmp %>%
    filter(PTM.CollapseKey %in% phospho_list_70)

  first_class_100 <- first_class_tmp %>%
    filter(PTM.CollapseKey %in% phospho_list_100)


  ### Finally, recreate PTM.CollapseKey with "GeneName_ResiduePosition_Multiplicity"
  # Split by semicolon and keep only the first element
  clean_ptm_dataframe <- function(df) {
    # Create temporary cleaned gene column
    df$PG.Gene_cleaned <- sapply(strsplit(df$PG.Genes, ";"), `[`, 1)

    # Build PTM.CollapseKey step by step
    df$PTM.CollapseKey <- paste0(
      df$PG.Gene_cleaned, "_",
      df$PTM.SiteAA,
      df$PTM.SiteLocation,
      "_M", df$PTM.Multiplicity
    )

    # Remove temporary column
    df <- df[, setdiff(names(df), "PG.Gene_cleaned")]

    # Remove duplicate PTM.CollapseKey rows
    df <- df[!duplicated(df$PTM.CollapseKey), ]

    # 5️⃣ Reset row names
    rownames(df) <- NULL

    return(df)
  }


  first_class_unfiltered <- clean_ptm_dataframe(first_class_unfiltered)
  first_class_50 <- clean_ptm_dataframe(first_class_50)
  first_class_70 <- clean_ptm_dataframe(first_class_70)
  first_class_100 <- clean_ptm_dataframe(first_class_100)

  ######################################################################################################
  ################ SECOND PART: CREATE INTENSITY/ANNOTATION/MMETADATA DAlist Object ####################
  ######################################################################################################

  make_DAList <- function(phosphosite_df, metadata) {
    # Ensure PG.Genes is cleaned and unique key built
    phosphosite_df$PG.Gene_cleaned <- sapply(strsplit(phosphosite_df$PG.Genes, ";"), `[`, 1)
    phosphosite_df$PTM.CollapseKey <- paste0(
      phosphosite_df$PG.Gene_cleaned, "_",
      phosphosite_df$PTM.SiteAA,
      phosphosite_df$PTM.SiteLocation,
      "_M", phosphosite_df$PTM.Multiplicity
    )
    phosphosite_df <- phosphosite_df[!duplicated(phosphosite_df$PTM.CollapseKey), ]
    phosphosite_df <- phosphosite_df[, setdiff(names(phosphosite_df), "PG.Gene_cleaned")]
    rownames(phosphosite_df) <- NULL

    # Identify quantity columns
    quantity_cols <- grep("PTM\\.Quantity", names(phosphosite_df), value = TRUE)

    # Build intensity and annotation data
    intensity_data <- phosphosite_df %>%
      dplyr::select(all_of(c(quantity_cols, "PTM.CollapseKey"))) %>%
      tibble::column_to_rownames("PTM.CollapseKey")

    annotation_data <- phosphosite_df %>%
      dplyr::select(!all_of(quantity_cols)) %>%
      tibble::column_to_rownames("PTM.CollapseKey") %>%
      as.data.frame()

    # Add uniprot_id for DAList compatibility
    annotation_data$uniprot_id <- rownames(annotation_data)

    # Simplify column names using metadata
    colnames(intensity_data) <- metadata$Simplify_ID
    rownames(metadata) <- metadata$Simplify_ID

    # Rename site probability and PTM group columns
    siteprob_cols <- grep("PTM\\.SiteProbability", names(annotation_data))
    PTMgroup_cols <- grep("PTM\\.Group", names(annotation_data))

    if (length(siteprob_cols) == nrow(metadata)) {
      colnames(annotation_data)[siteprob_cols] <-
        paste0(metadata$Simplify_ID, "_PTM_SiteProbability")
    }

    if (length(PTMgroup_cols) == nrow(metadata)) {
      colnames(annotation_data)[PTMgroup_cols] <-
        paste0(metadata$Simplify_ID, "_PTM_Group")
    }

    # Construct DAList
    raw <- proteoDA::DAList(
      data = intensity_data,
      annotation = annotation_data,
      metadata = metadata,
      design = NULL,
      eBayes_fit = NULL,
      results = NULL,
      tags = NULL
    )

    return(raw)
  }


  raw_unfiltered <- make_DAList(first_class_unfiltered, metadata)
  raw_50pc       <- make_DAList(first_class_50, metadata)
  raw_70pc       <- make_DAList(first_class_70, metadata)
  raw_100pc      <- make_DAList(first_class_100, metadata)


  #################################################################################
  ############## THIRD STEP: FILTER FOR VALID VALUE ###############################
  #################################################################################


  phosphosite_unfilter <- raw_unfiltered

  phosphosite_50pc <- filter_proteins_by_proportion_onegrp(raw_50pc,
                                                           min_prop = 0.50,
                                                           grouping_column = "group")

  phosphosite_70pc <- filter_proteins_by_proportion_onegrp(raw_70pc,
                                                           min_prop = 0.70,
                                                           grouping_column = "group")

  phosphosite_100pc<- filter_proteins_by_proportion_onegrp(raw_100pc,
                                                           min_prop = 1,
                                                           grouping_column = "group")

  list(
    phosphosite_unfilter = phosphosite_unfilter,
    phosphosite_50pc = phosphosite_50pc,
    phosphosite_70pc = phosphosite_70pc,
    phosphosite_100pc = phosphosite_100pc
  )


}

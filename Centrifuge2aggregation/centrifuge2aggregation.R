process_centrifuge_files <- function(directory_path, ncbi_db_path, filepattern) {
  # Get all files ending with "pattern"
  files <- list.files(directory_path, pattern = filepattern, full.names = TRUE)
  
  for (file_path in files) {
    cat("Processing file:", basename(file_path), "\n")
    
    # Read centrifuge report
    centrifuge_report <- read_tsv(file_path)
    
    # Taxonomy search
    taxID_to_search <- centrifuge_report %>% select(taxID) %>%
      unique() %>% pull()
    taxonomy_res <- as.data.frame(getTaxonomy(taxID_to_search, ncbi_db_path))
    taxonomy_res <- taxonomy_res %>% 
      mutate(taxID = taxID_to_search) %>%
      relocate(taxID) %>% 
      mutate(superkingdom = ifelse(taxID == 1, 'root', superkingdom))
    
    parsed_centrifuge <- dplyr::full_join(centrifuge_report, taxonomy_res, by = 'taxID')
    
    # Group unclassified
    na_superkingdom <- parsed_centrifuge %>% 
      filter(is.na(superkingdom))
    
    # Write excluded taxa
    excluded_file <- file.path(directory_path, paste0(tools::file_path_sans_ext(basename(file_path)), "_excluded_taxa.tsv"))
    write_tsv(na_superkingdom, excluded_file)
    
    na_superkingdom <- na_superkingdom %>% pull(name)
    
    if (length(na_superkingdom) > 0) {
      warning(paste("For file", basename(file_path), "the following entries have NA in superkingdom and will be considered as unclassified:",
                    paste(na_superkingdom, collapse = ", ")))
    }
    
    # Sum numReads and numUniqueReads for NA superkingdom entries
    na_sums <- parsed_centrifuge %>%
      filter(is.na(superkingdom)) %>%
      summarise(
        additional_numReads = sum(numReads),
        additional_numUniqueReads = sum(numUniqueReads)
      )
    
    centrifuge_processed <- parsed_centrifuge %>%
      mutate(
        numReads = if_else(name == "root", 
                           numReads + na_sums$additional_numReads, 
                           numReads),
        numUniqueReads = if_else(name == "root", 
                                 numUniqueReads + na_sums$additional_numUniqueReads, 
                                 numUniqueReads)
      ) %>%
      filter(!is.na(superkingdom))
    
    # Define taxonomic levels and find lowest level
    tax_levels <- c("species", "genus", "family", "order", "class", "phylum", "superkingdom")
    
    find_lowest_level <- function(row) {
      for (level in rev(tax_levels)) {
        if (!is.na(row[[level]])) {
          return(level)
        }
      }
      return(NA)
    }
    
    df_result <- centrifuge_processed %>%
      rowwise() %>%
      mutate(aggregation_level = find_lowest_level(cur_data())) %>%
      ungroup() %>%
      mutate(aggregated_counts = NA_real_)
    
    # Aggregate counts
    for (level in tax_levels) {
      df_result <- df_result %>%
        group_by(.data[[level]]) %>%
        mutate(tmp_agg = sum(numReads)) %>%
        mutate(aggregated_counts = if_else(aggregation_level == level, tmp_agg, aggregated_counts)) %>%
        ungroup()
    }
    
    # Remove the temporary aggregation column
    df_result <- df_result %>%
      select(-tmp_agg)
    
    # Write parsed result
    output_file <- file.path(directory_path, paste0(tools::file_path_sans_ext(basename(file_path)), "_parsed_centrifuge.tsv"))
    write_tsv(df_result, output_file)
    
    cat("Parsed report saved as:", basename(output_file), "\n\n")
  }
}

# Usage example:
# Replace with your actual directory and database paths
# directory_path <- "path/to/your/directory"
# ncbi_db_path <- "path/to/your/TaxonomizR_db/NCBItax_25.07.24.sql"
# filepattern <- "_report.txt"
# process_centrifuge_files(directory_path, ncbi_db_path, filepattern)
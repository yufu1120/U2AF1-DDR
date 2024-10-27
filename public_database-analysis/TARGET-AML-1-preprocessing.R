
# this script aims to load the counts (FPKM or TPM) from the MDACC bioinformatic core:
# https://bioinformatics.mdanderson.org/public-datasets/
{
  
  library(tibble)
  library(dplyr)
  library(readr)
  library(openxlsx)
  library(ggplot2)
  
  
  dir.path<-r'(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\AML-public_database\Data_from_MDACC)'
  setwd(dir.path)
  list.files()
  
  
  # 1. process the gene expression 
  
  # FPKM
  {
  cohort<-"TARGET-MDACC-FPKM"

  
  file<-r'(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\AML-public_database\Data_from_MDACC\424b6ab7e1ad5d0d6a8b2edd94dffc2a-data(TARGET-AML-FPKM)\original\original_matrix_data.tsv)'
  count.raw<-read_tsv(file)
  colnames(count.raw)
  count.raw <- count.raw %>% separate("...1", into = c("symbol", "ENSG"), sep = "\\|")
  
  # Remove rows that contain "PAR_Y" in the ENSG column
  count_par_y <- count.raw %>%
    filter(grepl("PAR_Y", ENSG))
  
  count_filtered_y <- count.raw %>%
    filter(!grepl("PAR_Y", ENSG))
  
  # Step 1: Identify repeated gene symbols
  repeated_genes <- count_filtered_y %>%
    group_by(symbol) %>%
    filter(n() > 1) %>%
    distinct(symbol) %>%
    pull(symbol)
  
  # Step 2: Filter the dataframe to separate the repeated and non-repeated rows
  count_repeats <- count_filtered_y %>%
    filter(symbol %in% repeated_genes)
  
  count_non_repeats <- count_filtered_y %>%
    filter(!symbol %in% repeated_genes)
  
  # Step 3: Calculate the mean for the repeated rows
  count_avg_repeats <- count_repeats %>%
    group_by(symbol) %>%
    summarize(across(where(is.numeric), mean, na.rm = T))
  
  # Step 4: Combine the unique rows and the averaged rows
  count_final <- bind_rows(count_non_repeats, count_avg_repeats)
  
  # Print the final dataframe
  print(count_final)
  count <- count_final[,3:ncol(count_final)]
  rownames(count) <- count_final$symbol
  colnames(count)
  # Remove "aq-" from the column names
  colnames(count) <- gsub("aq-", "", colnames(count))
  count.FPKM<-count
  }
  
  # then process the TPM=========================================================
  {
    cohort<-"TARGET-MDACC-TPM"
    # load the gene expression and metadata
    library(readr)
    file<-r'(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\AML-public_database\Data_from_MDACC\2f2120e54a9da9fe233ca231b4d5eac8-data(TARGET-AML-TPM)\original\original_matrix_data.tsv)'
    count.raw<-read_tsv(file)
    colnames(count.raw)
    count.raw$...1
    count.raw <- count.raw %>% separate("...1", into = c("symbol", "ENSG"), sep = "\\|")
    
    library(dplyr)
    # Remove rows that contain "PAR_Y" in the ENSG column
    count_par_y <- count.raw %>%
      filter(grepl("PAR_Y", ENSG))
    
    count_filtered_y <- count.raw %>%
      filter(!grepl("PAR_Y", ENSG))
    
    # Step 1: Identify repeated gene symbols
    repeated_genes <- count_filtered_y %>%
      group_by(symbol) %>%
      filter(n() > 1) %>%
      distinct(symbol) %>%
      pull(symbol)
    
    # Step 2: Filter the dataframe to separate the repeated and non-repeated rows
    count_repeats <- count_filtered_y %>%
      filter(symbol %in% repeated_genes)
    
    count_non_repeats <- count_filtered_y %>%
      filter(!symbol %in% repeated_genes)
    
    # Step 3: Calculate the mean for the repeated rows
    count_avg_repeats <- count_repeats %>%
      group_by(symbol) %>%
      summarize(across(where(is.numeric), mean, na.rm = T))
    
    # Step 4: Combine the unique rows and the averaged rows
    count_final <- bind_rows(count_non_repeats, count_avg_repeats)
    
    # Print the final dataframe
    print(count_final)
    count <- count_final[,3:ncol(count_final)]
    rownames(count) <- count_final$symbol
    colnames(count)
    # Remove "aq-" from the column names
    colnames(count) <- gsub("aq-", "", colnames(count))
    count.TPM<-count
  }
  
  # save the count here:
  cohort<-"TARGET-MDACC-FPKM"
  file.name<-paste(Sys.Date(),cohort,"rawdata-count.RDS", sep = "-")
  print(file.name)
  # saveRDS(count.FPKM, file = file.name)
  
  cohort<-"TARGET-MDACC-TPM"
  file.name<-paste(Sys.Date(),cohort,"rawdata-count.RDS", sep = "-")
  print(file.name)
  # saveRDS(count.TPM, file = file.name)
  
  # 2. the metadata has been processed somewhere else

  

  # load the FPKM or TPM from here:
  list.files()
  count<-readRDS("2024-09-29-TARGET-MDACC-FPKM-rawdata-count.RDS")
  count<-readRDS("2024-09-29-TARGET-MDACC-TPM-rawdata-count.RDS")
  clinical.metadata<-readRDS("clinical.metadata.merge.RDS")
  sample.metadata<-readRDS("sample.metadata.RDS")

}




# this script aims to load the counts (FPKM or TPM) from the MDACC bioinformatic core:
# https://bioinformatics.mdanderson.org/public-datasets/
{
  
  library(tibble)
  library(tidyr)
  library(readr)
  library(dplyr)
  library(openxlsx)
  library(ggplot2)
  
  

  dir.path<-r'(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\AML-public_database\Data_from_MDACC)'
  setwd(dir.path)
  list.files()
  
 
  # 1. process the gene expression
  
  {
  cohort<-"BEAT_AML1.0-MDACC-FPKM"
  file<-r'(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\AML-public_database\Data_from_MDACC\91d2a9afe9ca018c11896fd5caff6d42-data(BEAT-AML1.0-FPKM)\original\original_matrix_data.tsv)'
  count.raw<-read_tsv(file)
  colnames(count.raw)
  count.raw$...1
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
  cohort<-"BEAT_AML1.0-MDACC-TPM"
  
  file<-r'(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\AML-public_database\Data_from_MDACC\33592bfbb04a567137a02d352eabb6c3-data(BEAT-AML1_0-TPM)\original\original_matrix_data.tsv)'
  count.raw<-read_tsv(file)
  colnames(count.raw)
  count.raw$...1
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
  count.TPM<-count
  }
  
  # save the count here:
  cohort<-"BEAT_AML1.0-MDACC-FPKM"
  file.name<-paste(Sys.Date(),cohort,"rawdata-count.RDS", sep = "-")
  print(file.name)
  # saveRDS(count.FPKM, file = file.name)
  
  cohort<-"BEAT_AML1.0-MDACC-TPM"
  file.name<-paste(Sys.Date(),cohort,"rawdata-count.RDS", sep = "-")
  print(file.name)
  # saveRDS(count.TPM, file = file.name)
  
  
  
  
  
  # 2. process the metadata
  
  {
    metadata <- read.xlsx(xlsxFile = "beataml_wv1to4_clinical 1.xlsx",sheet = 1,
                          skipEmptyRows = FALSE,colNames=F)
    colnames(metadata)<-metadata[1,]
    metadata<-metadata[c(2:nrow(metadata)),]

    # saveRDS(metadata,"beataml_clinical.RDS")
  }
  
  
  # load the files from here:
  list.files()
  count<-readRDS("2024-09-29-BEAT_AML1.0-MDACC-FPKM-rawdata-count.RDS")
  count<-readRDS("2024-09-29-BEAT_AML1.0-MDACC-TPM-rawdata-count.RDS")
  metadata<-readRDS("beataml_clinical.RDS")

  
}



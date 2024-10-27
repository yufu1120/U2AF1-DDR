
# code 0: extract the output from rMATS and merge them into a summary table

# setwd to the rMATS output folder
setwd()
list.files()

# package
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(gplots)
library(ggvenn)
library(writexl)


#read all files
project_name<-"20230906_HDAC8-134"


AS_type_all<-c("SE","RI","MXE","A3SS","A5SS")
dir_all<-c("1","2")
dir_number<-length(dir_all)
cohort_all<-c("Ctrl-IR vs Ctrl-NIR","CM-IR vs CM-NIR")

# to read all files from rMATS

for (j in 1:dir_number) {
  dir<-dir_all[j]
  cohort<-cohort_all[j]
  
  for ( i in 1:length(AS_type_all)){
    AS_type<-AS_type_all[i]
    table_1<-read.table(paste0(".\\",project_name,"\\", dir,"\\final_output\\", AS_type,".MATS.JCEC.txt"), sep = "\t", header = T)
    table_1$cohort<-cohort
    table_1$AS_type<-AS_type
    
    if (AS_type=="SE"){
      table_1$feature_label<-paste0(table_1$chr,":",table_1$upstreamEE,",",table_1$exonStart_0base,"-",table_1$exonEnd,",",table_1$downstreamES)
      table_2<-table_1 %>%
        select(ID,GeneID,cohort,AS_type,feature_label,strand,PValue,FDR,IncLevel1,IncLevel2,IncLevelDifference)
    }
    
    if (AS_type=="RI"){
      table_1$feature_label<-paste0(table_1$chr,":",table_1$upstreamEE,",",table_1$riExonStart_0base,"-",table_1$riExonEnd,",",table_1$downstreamES)
      table_2<-table_1 %>%
        select(ID,GeneID,cohort,AS_type,feature_label,strand,PValue,FDR,IncLevel1,IncLevel2,IncLevelDifference)
    }
    
    if (AS_type=="MXE"){
      table_1$feature_label<-paste0(table_1$chr,":",table_1$upstreamEE,",",table_1$X1stExonStart_0base,"-",table_1$X1stExonEnd,",",table_1$X2ndExonStart_0base,",",table_1$X2ndExonEnd,",",table_1$downstreamES)
      table_2<-table_1 %>%
        select(ID,GeneID,cohort,AS_type,feature_label,strand,PValue,FDR,IncLevel1,IncLevel2,IncLevelDifference)
    }
    
    if (AS_type=="A3SS"){
      table_1$feature_label<-paste0(table_1$chr,":",table_1$longExonStart_0base,",",table_1$longExonEnd,",",table_1$shortES,",",table_1$shortEE,",",table_1$flankingES,",",table_1$flankingEE)
      table_2<-table_1 %>%
        select(ID,GeneID,cohort,AS_type,feature_label,strand,PValue,FDR,IncLevel1,IncLevel2,IncLevelDifference)
      
    }
    
    if (AS_type=="A5SS"){
      table_1$feature_label<-paste0(table_1$chr,":",table_1$longExonStart_0base,",",table_1$longExonEnd,",",table_1$shortES,",",table_1$shortEE,",",table_1$flankingES,",",table_1$flankingEE)
      table_2<-table_1 %>%
        select(ID,GeneID,cohort,AS_type,feature_label,strand,PValue,FDR,IncLevel1,IncLevel2,IncLevelDifference)
    }
    
    assign(paste(AS_type),table_2)
    
  }

# merge tables
summary_table<-as.data.frame(matrix(nrow =0, ncol = ncol(SE)))
colnames(summary_table)<-colnames(SE)

summary_table<-bind_rows(SE,RI, MXE, A5SS, A3SS) %>%
  separate(IncLevel1, c('IncLevel1.1', 'IncLevel1.2','IncLevel1.3'), sep = ",") %>%
  separate(IncLevel2, c('IncLevel2.1', 'IncLevel2.2','IncLevel2.3'), sep = ",") %>%
  mutate_at(c(9:14), as.numeric) %>%
  rowwise() %>%
  mutate(mean1 = mean(c_across(IncLevel1.1:IncLevel1.3),na.rm = T)) %>%
  mutate(mean2 = mean(c_across(IncLevel2.1:IncLevel2.3),na.rm = T)) %>%
  mutate(feature_label_short = paste0(AS_type,"-",ID,"-",GeneID)) %>%
  relocate(feature_label_short, .after = feature_label)

assign(cohort,summary_table)
setwd(paste0(".\\",project_name))
saveRDS(summary_table,paste0(Sys.Date(),'-splicing_summary-',cohort,'.RDS'))
write_xlsx(summary_table,paste0(Sys.Date(),'-splicing_summary-',cohort,'.xlsx'))
setwd("C:\\Users\\yufu\\OneDrive - City of Hope National Medical Center\\rotation\\1_Dr. Kuo\\U2AF1\\RNAseq\\splicing\\rMATS\\")

}

summary_table<-bind_rows(get(cohort_all[1]), get(cohort_all[2]))
saveRDS(summary_table, paste0(Sys.Date(),'-splicing_summary-all','.RDS'))
write_xlsx(summary_table,paste0(Sys.Date(),'-splicing_summary-all','.xlsx'))

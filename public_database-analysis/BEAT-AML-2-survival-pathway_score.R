
# this script meant to use the harmonized FPKM counts from the website (https://biodev.github.io/BeatAML2/),
# to compare the gene expression between AML and healthy, as well as calculating the pathway score by averaging the expression
{
  
  library(tibble)
  library(tidyr)
  library(dplyr)
  library(readr)
  library(openxlsx)
  library(ggplot2)
  
  
  cohort<-"BEAT_AML1.0-MDACC-FPKM"
  # cohort<-"BEAT_AML1.0-MDACC-TPM"
  dir.path<-r'(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\AML-public_database\Data_from_MDACC)'
  setwd(dir.path)
  list.files()
  
 

  # load the files from here:
  count<-readRDS("2024-09-29-BEAT_AML1.0-MDACC-FPKM-rawdata-count.RDS")
  # count<-readRDS("2024-09-29-BEAT_AML1.0-MDACC-TPM-rawdata-count.RDS")
  metadata<-readRDS("beataml_clinical.RDS")
  
  # read gene set file
  xlsxFile.path<-r"(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\AML-public_database\fimmu.2021.806324-Table-2.xlsx)"
  signature.file <- openxlsx::read.xlsx(xlsxFile = xlsxFile.path ,sheet = 1, skipEmptyRows = FALSE, colNames=T)
  
}


{
  #  include the expression of gene of interest 

# define the gene of interest
signature.genes <- signature.file[9]
# gene.symbols1 <- signature.genes[1:nrow(signature.genes),]
# gene.symbols.name1<-"DDR_score"
gene.symbols1 <- ""
gene.symbols.name1<-""

###### calculate the pathway score
intersect.genes1 <- intersect(gene.symbols1, row.names(count))
genes.notappear1 <- setdiff(gene.symbols1,intersect.genes1) # genes not included in gene expression data
extrgenes.matr1 <- count[intersect.genes1,]
gene_expression1 <- as.matrix(colMeans(extrgenes.matr1,na.rm=T)) %>%
  as.data.frame()
colnames(gene_expression1) <- gene.symbols.name1

# and define the gene of interest
gene.symbols2 <-"HDAC8"
gene.symbols.name2<-gene.symbols2
###### calculate the pathway score
intersect.genes2 <- intersect(gene.symbols2, row.names(count))
genes.notappear2 <- setdiff(gene.symbols2,intersect.genes2) # genes not included in gene expression data
extrgenes.matr2 <- count[intersect.genes2,]
gene_expression2 <- as.matrix(colMeans(extrgenes.matr2,na.rm=T)) %>%
  as.data.frame()
colnames(gene_expression2) <- gene.symbols.name2
gene_expression <-left_join(rownames_to_column(gene_expression1), rownames_to_column(gene_expression2), by=c("rowname" = "rowname") )

GOI <-c("HDAC8","IL1RL1")
intersect.genes2 <- intersect(GOI, row.names(count))
genes.notappear2 <- setdiff(GOI,intersect.genes2) # genes not included in gene expression data
gene_expression<-count[rownames(count) %in% GOI, ] %>% 
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "patient_id")

colnames(gene_expression)<-c("patient_id", GOI)

} 

# an function to calculate the pathway score
cal_pathway_score<-function(signature.file, count){
  empty.matx<-data.frame()
  for (i in 1:ncol(signature.file)){
    gene.symbols1 <- signature.file[i] %>%
      na.exclude() %>% t()
    gene.symbols.name1 <- colnames(signature.file)[i]
    ###### calculate the pathway score
    intersect.genes1 <- intersect(gene.symbols1, row.names(count))
    # genes.notappear1 <- setdiff(gene.symbols1,intersect.genes1) # genes not included in gene expression data
    extrgenes.matr1 <- count[intersect.genes1,]
    gene_expression1 <- as.matrix(colMeans(extrgenes.matr1,na.rm=T)) %>%
      as.data.frame()
    colnames(gene_expression1) <- gene.symbols.name1
    
    if(ncol(empty.matx) == 0){
      empty.matx = gene_expression1
      
    } else{
      empty.matx = cbind(empty.matx, gene_expression1)
    }
    rownames(empty.matx) = colnames(count)
  }
  empty.matx<-empty.matx %>%
    rownames_to_column()
  return(empty.matx)
}

score_summary<-cal_pathway_score(signature.file, count)

# merge two gene score and metadata
{
  df<-gene_expression %>%
  left_join(., score_summary,  by=c("patient_id" = "rowname")) %>%
  left_join(., metadata, by=c("patient_id" = "dbgap_rnaseq_sample")) %>%                      #left 707
  dplyr::filter(used_manuscript_analyses=="yes") %>%
  dplyr::filter(diseaseStageAtSpecimenCollection %in% c("Initial Diagnosis", "Healthy"))
  
  # filter any specific subtypes if necessary:
  # dplyr::filter(consensusAMLFusions=="CBFB-MYH11")
  table(df$consensusAMLFusions, useNA="always")
  
  df<-df %>%
    dplyr::mutate(Gene_Fusion_special = dplyr::case_when(
      diseaseStageAtSpecimenCollection == "Healthy"  ~ "Healthy",
      diseaseStageAtSpecimenCollection != "Healthy"  ~ consensusAMLFusions)) %>%
    relocate(Gene_Fusion_special, .before = 2) %>%
    relocate(diseaseStageAtSpecimenCollection, .after = 2) %>%
    relocate(specimenType, .after = 3)
  
  df$Gene_Fusion_special[is.na(df$Gene_Fusion_special)]<-"Unknown"
  table(df$Gene_Fusion_special, useNA="always")
  
  # write table

  file.name<-paste(Sys.Date(),cohort,"rawdata-DDR_score(with_healthy).xlsx", sep = "-")
  print(file.name)
  # write.xlsx(df,file =file.name)
  
}

# check DDR in different variable
png(paste(Sys.Date(),cohort,"box plot-fusion-DDR",".png",sep="-"), res=300, units="in", height=6, width=6)
colnames(df)
unique(df$specimenType)
df%>%
  # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
  ggplot(., aes(Gene_Fusion_special, DDR_human_GO00006974, fill=Gene_Fusion_special, color=specimenType)) +
  geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
graphics.off()  
# ggpubr::stat_compare_means(method = "t.test")
# ggpubr::stat_compare_means()
compare_means(HR ~ consensusAMLFusions, data = df)%>%
  print(n=21)

# check HDAC8 in different variable
png(paste(Sys.Date(),cohort,"box plot-fusion-HDAC8",".png",sep="-"), res=300, units="in", height=6, width=6)
colnames(df)
unique(df$specimenType)
df%>%
  # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
  ggplot(., aes(Gene_Fusion_special, HDAC8, fill=Gene_Fusion_special, color=specimenType)) +
  geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
graphics.off()  
# ggpubr::stat_compare_means(method = "t.test")
# ggpubr::stat_compare_means()
compare_means(HR ~ consensusAMLFusions, data = df)%>%
  print(n=21)





{ # further filter samples with survival information to make survival curve
  df <- df %>% 
    dplyr::mutate(status = dplyr::case_when(
    vitalStatus == "Dead"~ 1,
    vitalStatus == "Unknown"~ NA,
    vitalStatus == "Alive"~ 0)) %>%
    dplyr::filter(!is.na(status))%>%
    dplyr::filter(!is.na(overallSurvival))
  df$overallSurvival<-as.numeric(df$overallSurvival)
  
}

# determine the survival cutpoint:
# using the maximally selected rank statistics from the 'maxstat' R package
# {
#   # define cut for gene 1
#   res.cut1<-survminer::surv_cutpoint(
#     data=df,
#     time = "overallSurvival",
#     event = "status",
#     variables = gene.symbols.name1,
#     minprop = 0.1,
#     progressbar = TRUE)
#   
#   plot(res.cut1, gene.symbols.name1, palette = "png")
#   summary(res.cut1)
#   res.cat1 <- surv_categorize(res.cut1)
#   fit1 <- survfit(Surv(overallSurvival, status) ~ get(gene.symbols.name1), data = res.cat1)
#   
#   df<-df %>%
#     mutate(gene.symbols1_cat = res.cat1[,gene.symbols.name1]) %>%
#     relocate(gene.symbols1_cat, .after = gene.symbols.name1)
  
  # define cut for gene 2
  res.cut2<-survminer::surv_cutpoint(
    data=df,
    time = "overallSurvival",
    event = "status",
    variables = gene.symbols.name2,
    minprop = 0.1,
    progressbar = TRUE)
  
  file.name<-paste(Sys.Date(),cohort,"survival_cut",gene.symbols.name2,".png", sep = "-")
  png(file.name, res=300, units="in", height=6, width=5)
  
  print(file.name)
  plot(res.cut2, gene.symbols.name2, palette = "png")
  graphics.off()
  
  summary(res.cut2)
  res.cat2 <- surv_categorize(res.cut2)
  fit2 <- survfit(Surv(overallSurvival, status) ~ get(gene.symbols.name2), data = res.cat2)
  
  df<-df %>%
    mutate(gene.symbols2_cat = res.cat2[,gene.symbols.name2]) %>%
    relocate(gene.symbols2_cat, .after = gene.symbols.name2)
  
  # write table
  df.write<-df %>%
    relocate(overallSurvival, .after = gene.symbols2_cat) %>%
    relocate(status, .after = overallSurvival)
  file.name<-paste(Sys.Date(),cohort,"survival_rawdata",gene.symbols.name1,"GO",".xlsx", sep = "-")
  print(file.name)
  write.xlsx(df.write,file =file.name)
} 
  

# survfit(formula = lungsurvival ~ sex + ph.ecog, data=lung)
  {
    png(file.name, res=300, units="in", height=6, width=5)
    file.name<-paste(Sys.Date(),cohort,"survival_plot",gene.symbols.name1,".png", sep = "-")
    ggsurvplot(
      fit1,
      data = res.cat1,
      title  = paste0(gene.symbols.name1, " (BEAT-AML Waves-all)"),
      break.time.by = 500,
      pval = TRUE,                # Show p-value
      conf.int = F,            # Show confidence intervals
      risk.table = TRUE,          # Show risk table
      legend.title = "Gene Expression",
      # legend.labs = c(1:8),
      legend.labs = c("High","Low"),
      xlab = "Days",
      ylab = "Survival Probability"
    )
    graphics.off()
  }

{
  file.name<-paste(Sys.Date(),cohort,"survival_plot",gene.symbols.name2,".png", sep = "-")
  png(file.name, res=300, units="in", height=6, width=5)
  ggsurvplot(
    fit2,
    data = res.cat2,
    title  = paste0(gene.symbols.name2, " (BEAT-AML Waves-all)"),
    break.time.by = 500,
    pval = TRUE,                # Show p-value
    conf.int = F,            # Show confidence intervals
    risk.table = TRUE,          # Show risk table
    legend.title = "Gene Expression",
    # legend.labs = c(1:8),
    legend.labs = c("High","Low"),
    xlab = "Days",
    ylab = "Survival Probability"
  )
  graphics.off()
}

df$consensusAMLFusions

# bar plot
png(paste(Sys.Date(),cohort,"box plot-DDR-HDAC8",".png",sep="-"), res=300, units="in", height=6, width=6)
df%>%
  # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
  ggplot(., aes(gene.symbols1_cat, HDAC8, fill=gene.symbols1_cat, color=gene.symbols1_cat)) +
  geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))
graphics.off()  
# ggpubr::stat_compare_means(method = "t.test")
# ggpubr::stat_compare_means()
compare_means(HDAC8 ~ gene.symbols1_cat, data = df)
group.by = "gene_expression_cat")%>%
    print(n=25)
  
  # check variables in HDAC8 high vs low
  png(paste(Sys.Date(),cohort,"box plot-HDAC8-DDR",".png",sep="-"), res=300, units="in", height=6, width=6)
  df%>%
    # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
    ggplot(., aes(gene.symbols2_cat, DDR, fill=gene.symbols2_cat, color=gene.symbols2_cat)) +
    geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))+
    xlab("HDAC8")+
    labs(color="HDAC8", fill="HDAC8")
  graphics.off()  
  # ggpubr::stat_compare_means(method = "t.test")
  # ggpubr::stat_compare_means()
  compare_means(DDR ~ gene.symbols2_cat, data = df)
  group.by = "gene_expression_cat")%>%
    print(n=25)
  
  png(paste(Sys.Date(),cohort,"box plot-HDAC8-HR",".png",sep="-"), res=300, units="in", height=6, width=6)
  df%>%
    # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
    ggplot(., aes(gene.symbols2_cat, HR, fill=gene.symbols2_cat, color=gene.symbols2_cat)) +
    geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))+
    xlab("HDAC8")+
    labs(color="HDAC8", fill="HDAC8")
  graphics.off()  
  # ggpubr::stat_compare_means(method = "t.test")
  # ggpubr::stat_compare_means()
  compare_means(HR ~ gene.symbols2_cat, data = df)
  group.by = "gene_expression_cat")%>%
    print(n=25)
  
  png(paste(Sys.Date(),cohort,"box plot-HDAC8-NHEJ",".png",sep="-"), res=300, units="in", height=6, width=6)
  df%>%
    # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
    ggplot(., aes(gene.symbols2_cat, NHEJ, fill=gene.symbols2_cat, color=gene.symbols2_cat)) +
    geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))+
    xlab("HDAC8")+
    labs(color="HDAC8", fill="HDAC8")
  graphics.off()  
  # ggpubr::stat_compare_means(method = "t.test")
  # ggpubr::stat_compare_means()
  compare_means(NHEJ ~ gene.symbols2_cat, data = df)
  group.by = "gene_expression_cat")%>%
    print(n=25)
  
  # check DDR in different variable
  png(paste(Sys.Date(),cohort,"box plot-fusion-DDR",".png",sep="-"), res=300, units="in", height=6, width=6)
  df%>%
    # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
    ggplot(., aes(consensusAMLFusions, DDR, fill=consensusAMLFusions, color=consensusAMLFusions)) +
    geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  graphics.off()  
  # ggpubr::stat_compare_means(method = "t.test")
  # ggpubr::stat_compare_means()
  compare_means(HR ~ consensusAMLFusions, data = df)%>%
    print(n=21)
  
  
  # heatmap
  library(ComplexHeatmap)
  # scale (scale() applies by columns)
  df_order <-df %>%
    arrange(gene.symbols2_cat)
  mat<-df_order[,6:14] %>%
    scale()
  
  row_ha = rowAnnotation(HDAC8 = df_order$gene.symbols2_cat,
                         col = list(HDAC8 = c("high" = "#F8766D", "low" = "#00BFC4")))
  png(paste(Sys.Date(),cohort,"heatmap-HDAC8-DDR",".png",sep="-"), res=300, units="in", height=6, width=6)
  ComplexHeatmap::Heatmap(mat, right_annotation = row_ha, name = "pathway\nscore", width = unit(6, "cm"), column_names_rot = 45)
  graphics.off()
  
  ########### or Classify patients into high and low gene expression groups
  median_expression <- median(df$gene_expression)
  
  # define groups
  df$expression_group <- ifelse(df$gene_expression > median_expression  , "High", "Low")
  df$expression_group<-factor(df$expression_group, levels = c("High", "Low"))
  surv_object <- Surv(time = df$overallSurvival, event = df$status)
  # Fit the survival model
  # fit <- survfit(surv_object ~ expression_group+specimenType, data = df)
  fit <- survfit(surv_object ~ expression_group, data = df)
  
  # Plot the survival curve
  if(number_bins==2){
    ggsurvplot(
      fit,
      data = df,
      # title  = gene.symbols,
      break.time.by = 500,
      pval = TRUE,                # Show p-value
      # conf.int = TRUE,            # Show confidence intervals
      risk.table = TRUE,          # Show risk table
      legend.title = "Gene Expression",
      # legend.labs = c("High", "Low"),
      # legend.labs = c("WT", "Mutant"),
      xlab = "Days",
      ylab = "Survival Probability"
    )
  }
  ###########
  
  
  {
    df2<-gene_expression %>%
      as.data.frame() %>%
      dplyr::rename(patient_id=rowname)%>%
      left_join(., score_summary, , by=c("patient_id" = "rowname")) %>%
      left_join(., metadata, by=c("patient_id" = "dbgap_rnaseq_sample")) %>%                      #left 707
      dplyr::mutate(status = dplyr::case_when(
        vitalStatus == "Dead"~ 1,
        vitalStatus == "Unknown"~ NA,
        vitalStatus == "Alive"~ 0)) %>%
      dplyr::filter(!is.na(status))%>%
      dplyr::filter(!is.na(overallSurvival))%>%
      dplyr::filter(used_manuscript_analyses=="yes")
      # dplyr::filter(diseaseStageAtSpecimenCollection=="Initial Diagnosis") 
    df2$overallSurvival<-as.numeric(df2$overallSurvival)
  }
  
  png(paste(Sys.Date(),cohort,"box plot-diseaseStageAtSpecimenCollection-HDAC8",".png",sep="-"), res=300, units="in", height=6, width=6)
  df2%>%
    # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
    ggplot(., aes(diseaseStageAtSpecimenCollection, HDAC8, fill=diseaseStageAtSpecimenCollection, color=diseaseStageAtSpecimenCollection)) +
    geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))
  graphics.off()  
  # ggpubr::stat_compare_means(method = "t.test")
  # ggpubr::stat_compare_means()
  compare_means(HDAC8 ~ diseaseStageAtSpecimenCollection, data = df2)
  group.by = "gene_expression_cat")%>%
    print(n=25)
  graphics.off()
  
  png(paste(Sys.Date(),cohort,"box plot-diseaseStageAtSpecimenCollection-DDR",".png",sep="-"), res=300, units="in", height=6, width=6)
  df2%>%
    # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
    ggplot(., aes(diseaseStageAtSpecimenCollection, DDR, fill=diseaseStageAtSpecimenCollection, color=diseaseStageAtSpecimenCollection)) +
    geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))
  graphics.off()  
  # ggpubr::stat_compare_means(method = "t.test")
  # ggpubr::stat_compare_means()
  compare_means(DDR ~ diseaseStageAtSpecimenCollection, data = df2)
  graphics.off()
  group.by = "gene_expression_cat")%>%
    print(n=25)

  length(unique(df$patient_id))
duplicated(df$patient_id)


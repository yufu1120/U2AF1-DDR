

{
library(IOBR)
library(tibble)
library(dplyr)
library(openxlsx)
library(haemdata)
library(SummarizedExperiment)

# load the gene expression and metadata
  
cohort<-"CM-2016-2018-BM"
dir.path<-r'(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\AML-public_database\CM-2016-TPM)'
setwd(dir.path)
list.files()

use_pinboard("onedrive")
haemdata::get_pin_list()

all_mice.metadata <- get_pin("metadata_mmu.csv")
metadata <-all_mice.metadata %>%
  filter(cohort == "AML.mRNA.2020") %>%
  filter(assay == "mRNA") %>%
  filter(tissue == "BM") %>%
  distinct(., library_id, .keep_all = T)

file.name<-paste(Sys.Date(),cohort,"metadata",".xlsx", sep = "-")
print(file.name)
write.xlsx(metadata,file =file.name)


all_mice <- get_pin("mmu_mrna_all_mice_GENCODEm28_HLT_qc.rds")
all_mice.count<-assays(all_mice)$counts
total.exprmatr<-all_mice.count[, colnames(all_mice.count) %in% metadata$library_id]
row.names(total.exprmatr)<- sub("\\..*", "", row.names(total.exprmatr))


xlsxFile.path<-r"(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\AML-public_database\mouse_DDR-pathway-gene-list.xlsx)"

# convert raw counts to TPM
tpm<-as.data.frame(IOBR::count2tpm(countMat = total.exprmatr, source = "local", idType = "Ensembl", org="mmus"))
}


{
# read gene set file
signature.file <- openxlsx::read.xlsx(xlsxFile = xlsxFile.path ,sheet = 1, skipEmptyRows = FALSE, colNames=T)
colnames(signature.file)
# covert Human to Mouse:
# {
#   convert_table<-read.csv("human_mouse_1to1_orthologs.csv")
# signature.file.converted<-c()
# for (i in 1:ncol(signature.file)){
#   signature.file.converted[[i]] <- convert_table$mouse[match(signature.file[[i]], convert_table$human)]
# 
# }
# signature.file.converted<-as.data.frame(signature.file.converted)
# colnames(signature.file.converted)<-colnames(signature.file)
# signature.file<-signature.file.converted
# }

# define the gene set
  signature.genes <- signature.file[1]
# signature.genes <- signature.file[9]
gene.symbols1 <- signature.genes %>% t()
gene.symbols.name1<-"DDR_score"
# gene.symbols.name1<-"top1_score"
###### calculate the pathway score
intersect.genes1 <- intersect(gene.symbols1, row.names(tpm))
genes.notappear1 <- setdiff(gene.symbols1,intersect.genes1) # genes not included in gene expression data
extrgenes.matr1 <- tpm[intersect.genes1,]
gene_expression1 <- as.matrix(colMeans(extrgenes.matr1,na.rm=T)) %>%
  as.data.frame()
colnames(gene_expression1) <- gene.symbols.name1

# and define the gene of interest
gene.symbols2 <-"Il1rl1"
gene.symbols.name2<-gene.symbols2
###### calculate the pathway score
intersect.genes2 <- intersect(gene.symbols2, row.names(tpm))
genes.notappear2 <- setdiff(gene.symbols2,intersect.genes2) # genes not included in gene expression data
extrgenes.matr2 <- tpm[intersect.genes2,]
gene_expression2 <- as.matrix(colMeans(extrgenes.matr2,na.rm=T)) %>%
  as.data.frame()
colnames(gene_expression2) <- gene.symbols.name2


gene_expression <-left_join(rownames_to_column(gene_expression1), rownames_to_column(gene_expression2), by=c("rowname" = "rowname") )
} 

# an function to calculate the pathway score
cal_pathway_score<-function(signature.file, tpm){
  empty.matx<-data.frame()
  for (i in 1:ncol(signature.file)){
    gene.symbols1 <- signature.file[i] %>%
      na.exclude() %>% t()
    gene.symbols.name1 <- colnames(signature.file)[i]
    ###### calculate the pathway score
    intersect.genes1 <- intersect(gene.symbols1, row.names(tpm))
    # genes.notappear1 <- setdiff(gene.symbols1,intersect.genes1) # genes not included in gene expression data
    extrgenes.matr1 <- tpm[intersect.genes1,]
    gene_expression1 <- as.matrix(colMeans(extrgenes.matr1,na.rm=T)) %>%
      as.data.frame()
    colnames(gene_expression1) <- gene.symbols.name1
    
    if(ncol(empty.matx) == 0){
      empty.matx = gene_expression1
      
    } else{
      empty.matx = cbind(empty.matx, gene_expression1)
    }
    rownames(empty.matx) = colnames(tpm)
  }
  empty.matx<-empty.matx %>%
    rownames_to_column()
  return(empty.matx)
}


# merge two gene score and metadata
{
  score_summary<-cal_pathway_score(signature.file, tpm)
  metadata$library_id<-as.character(metadata$library_id)
  df<-gene_expression %>%
  as.data.frame() %>%
  dplyr::rename(COH_ID=rowname)%>%
  left_join(., score_summary, by=c("COH_ID" = "rowname")) %>%
  left_join(., metadata, by=c("COH_ID" = "library_id")) %>%
  dplyr::arrange(treatment)
    
  
  # df$time<-factor(df$time, levels = c("T0",  "T1",   "T2",   "T3",   "T4",   "T5",   "T6",   "T6p5", "T7",   "T8",   "T9",   "T10"  ))
  # df$phenotype<-factor(df$phenotype, levels = c("CTL", "CM"))
  
  file.name<-paste(Sys.Date(),cohort,"raw_data",gene.symbols.name1,".xlsx", sep = "-")
  print(file.name)
  write.xlsx(df,file =file.name)
  

  }

# heatmap
library(ComplexHeatmap)
# scale (scale() applies by columns)


df2<-df %>%
  arrange(time, desc(genotype))

mat<-df2[,4:12] %>%
  scale()
row_ha = rowAnnotation(phenotype = df2$phenotype,
                       time=df2$time,
                       col = list(phenotype = c("CTL" = "gray", "CM" = "#00BFC4")))
ComplexHeatmap::Heatmap(mat, right_annotation = row_ha, cluster_rows = FALSE)


png(paste(Sys.Date(),cohort,"heatmap-HDAC8-DDR",".png",sep="-"), res=300, units="in", height=6, width=6)
ComplexHeatmap::Heatmap(mat, right_annotation = row_ha, name = "pathway\nscore", width = unit(6, "cm"), column_names_rot = 45)
graphics.off()


png(paste(Sys.Date(),cohort,"box plot-top1-time",".png",sep="-"), res=300, units="in", height=6, width=6)
df%>%
  # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
  ggplot(., aes(time, top1, fill=time, color=phenotype)) +
  geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))
graphics.off()  
# ggpubr::stat_compare_means(method = "t.test")
# ggpubr::stat_compare_means()
compare_means(HDAC8 ~ gene.symbols1_cat, data = df) %>%
    print(n=25)

  png(paste(Sys.Date(),cohort,"box plot-Il1rl1-time",".png",sep="-"), res=300, units="in", height=6, width=6)
  df.filtered<-df%>%
    filter(mouse_id %in% c("2690","2702","2705","2708","2718","2719","2720","2731")) 
    # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
  df.filtered %>%
    ggplot(., aes(treatment, DDR_mouse_GO00006974, fill=treatment, color=treatment)) +
    geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))
  graphics.off()  
  # ggpubr::stat_compare_means(method = "t.test")
  # ggpubr::stat_compare_means()
  compare_means(DDR_mouse_GO00006974 ~ treatment, data = df.filtered)
  group.by = "gene_expression_cat")%>%
    print(n=25)
  
  png(paste(Sys.Date(),cohort,"box plot-DDR-time",".png",sep="-"), res=300, units="in", height=6, width=6)
  df%>%
    # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
    ggplot(., aes(sample_date, DDR_mouse_GO00006974, fill=treatment, color=treatment)) +
    geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))
  graphics.off()  
  # ggpubr::stat_compare_means(method = "t.test")
  # ggpubr::stat_compare_means()
  compare_means(Il1rl1 ~ gene.symbols1_cat, data = df)
  group.by = "gene_expression_cat")%>%
    print(n=25)
  
  png(paste(Sys.Date(),cohort,"box plot-DDR-all-time",".png",sep="-"), res=300, units="in", height=6, width=6)
  df%>%
    # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
    ggplot(., aes(time, DDR_all, fill=time, color=phenotype)) +
    geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))
  graphics.off()  
  # ggpubr::stat_compare_means(method = "t.test")
  # ggpubr::stat_compare_means()
  df$time_phenotype<-paste(df$time, df$phenotype, sep = "_")
  table<-compare_means(DDR_all ~ time_phenotype , data = df)
  file.name<-paste(Sys.Date(),cohort,"stac",gene.symbols.name1,"all", "time", "phenotype",".xlsx", sep = "-")
  print(file.name)
  write.xlsx(table,file =file.name)

  
  png(paste(Sys.Date(),cohort,"box plot-DDR-small-time",".png",sep="-"), res=300, units="in", height=6, width=6)
  df%>%
    # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
    ggplot(., aes(time, DDR_small, fill=time, color=phenotype)) +
    geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))
  graphics.off()  
  # ggpubr::stat_compare_means(method = "t.test")
  # ggpubr::stat_compare_means()
  table<-compare_means(DDR_all ~ time_phenotype , data = df)
  file.name<-paste(Sys.Date(),cohort,"stac",gene.symbols.name1,"small", "time", "phenotype",".xlsx", sep = "-")
  print(file.name)
  write.xlsx(table,file =file.name)
  
  png(paste(Sys.Date(),cohort,"box plot-HR-time",".png",sep="-"), res=300, units="in", height=6, width=6)
  df%>%
    # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
    ggplot(., aes(time, HR, fill=time, color=phenotype)) +
    geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))
  graphics.off()  
  # ggpubr::stat_compare_means(method = "t.test")
  # ggpubr::stat_compare_means()
  compare_means(Il1rl1 ~ gene.symbols1_cat, data = df)
  group.by = "gene_expression_cat")%>%
    print(n=25)
  
  png(paste(Sys.Date(),cohort,"box plot-NHEJ-time",".png",sep="-"), res=300, units="in", height=6, width=6)
  df%>%
    # dplyr::filter(Risk_group%in%c("High Risk","Standard Risk","Low Risk"))%>%
    ggplot(., aes(time, NHEJ, fill=time, color=phenotype)) +
    geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))
  graphics.off()  
  # ggpubr::stat_compare_means(method = "t.test")
  # ggpubr::stat_compare_means()
  compare_means(Il1rl1 ~ gene.symbols1_cat, data = df)
  group.by = "gene_expression_cat")%>%
    print(n=25)





# determine the survival cutpoint:
# using the maximally selected rank statistics from the 'maxstat' R package
{
  # define cut for gene 1
  res.cut1<-survminer::surv_cutpoint(
    data=df,
    time = "overallSurvival",
    event = "status",
    variables = gene.symbols.name1,
    minprop = 0.1,
    progressbar = TRUE)
  
  plot(res.cut1, gene.symbols.name1, palette = "png")
  summary(res.cut1)
  res.cat1 <- surv_categorize(res.cut1)
  fit1 <- survfit(Surv(overallSurvival, status) ~ get(gene.symbols.name1), data = res.cat1)
  
  df<-df %>%
    mutate(gene.symbols1_cat = res.cat1[,gene.symbols.name1]) %>%
    relocate(gene.symbols1_cat, .after = gene.symbols.name1)
  
  # define cut for gene 2
  res.cut2<-survminer::surv_cutpoint(
    data=df,
    time = "overallSurvival",
    event = "status",
    variables = gene.symbols.name2,
    minprop = 0.1,
    progressbar = TRUE)
  
  plot(res.cut2, gene.symbols.name2, palette = "png")
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
  file.name<-paste(Sys.Date(),cohort,"survival_rawdata",gene.symbols.name1,gene.symbols.name2,".xlsx", sep = "-")
  print(file.name)
  # write.xlsx(df.write,file =file.name)
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
  png(file.name, res=300, units="in", height=6, width=5)
  file.name<-paste(Sys.Date(),cohort,"survival_plot",gene.symbols.name2,".png", sep = "-")
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


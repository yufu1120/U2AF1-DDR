analysis_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
analysis_output_dir <- file.path(analysis_dir, "output")
dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)

{
  # package
  library(devtools)
  library(ggbiplot)
  library(tidyverse)
  library(edgeR)
  library(ggrepel)
  library(sva)
  set.seed(126)
  
  
  # load raw data:
  data_file_1 <- "counts/new_counts_ykuo_221114_mm10.txt"
  data <- read.table(data_file_1, header = T, sep = "\t")
  colnames(data)
  data<-data[,c(6,8:ncol(data))]
  row.names(data)<-data[,1]
  data<-data[,2:ncol(data)]
  colnames(data) <- gsub("X","",colnames(data))
  
  
  
  metadata<-read.table("./metadata/metadata.txt",header = T)
  row.names(metadata)<-metadata$ID
  
  
  all(colnames(data) %in% rownames(metadata))
  index<-match(rownames(metadata),colnames(data))
  data<-data[,index]
  
  all(colnames(data) ==rownames(metadata))
  
  raw.count.all<-data
  metadata.all<-metadata
  
  # # save the metadata and raw count
  saveRDS(raw.count.all,"./counts/20230425_raw.count.all.rds")
  write.csv(raw.count.all,"./counts/20230425_raw.count.all.csv")
  
  colnames(metadata.all)
  
  metadata.all$cohort<-as.factor(metadata.all$cohort)
  
  metadata.all<- metadata.all %>% 
    unite(group2, c("phenotype", "treatment"), remove=F)
  
  saveRDS(metadata.all, "./metadata/20230425_metadata.all.rds")

  
  # make raw count into cpm ================================================
  # code from https://github.com/drejom/smrnaseq/blob/master/bin/edgeR_miRBase.r
  
  raw.to.cpm.filter<-function(data, meta){
    # Remove genes with 0 reads in all samples
    row_sub = apply(data, 1, function(row) all(row ==0 ))
    # Only subset if at least one sample is remaining
    nr_keep <- sum(row_sub)
    if (nr_keep > 0){
      data<-data[!row_sub,]
    }
    
    #Also check for colSums > 0, otherwise DGEList will fail if samples have entirely colSum == 0
    drop_colsum_zero <- (colSums(data, na.rm=T) != 0) # T if colSum is not 0, F otherwise
    data <- data[, drop_colsum_zero] # all the non-zero columns
    
    # Normalization
    dataDGE<-edgeR::DGEList(counts=data,genes=rownames(data), group=meta[['group']])
    o <- order(rowSums(dataDGE$counts), decreasing=TRUE)
    dataDGE <- dataDGE[o,]
    
    # TMM
    dataNorm <- edgeR::calcNormFactors(dataDGE)
    
    # filter out low expression genes
    y<-dataNorm
    keep <- edgeR::filterByExpr(y, design=NULL)
    y <- y[keep,]
    dataNorm <- y
    
    # Print normalized read counts to file
    dataNorm_df<-as.data.frame(edgeR::cpm(dataNorm))
    
    aa<-list("DEG"=dataNorm, "df"=dataNorm_df)
    return(aa)
    
  }
  
  cpm.all<-raw.to.cpm.filter(raw.count.all, metadata.all)
  colnames(cpm.all$df)<-metadata.all$simple_ID
  
  # saveRDS(cpm.all,"20230425_cpm.count.all.rds")
  # write.csv(cpm.all$df,"20230425_cpm.count.all.csv")
  
  cpm.all<-readRDS("./counts/20230425_cpm.count.all.rds")
  metadata.all<-readRDS("./metadata/20230425_metadata.all.rds")
  
  # split cohort:
  # Retain the 14 CM/WT samples listed in metadata; count columns refer to the
  cm_count_columns <- c(1, 3:9, 11:16)
  CHW.IR.cpm<-cpm.all$df[,cm_count_columns]
  CHW.IR.meta<-metadata.all
  stopifnot(ncol(CHW.IR.cpm) == nrow(CHW.IR.meta))
  
}


# PCA ==========================================================================
PCA_plot<-function(a,meta,suffix){
  # perfrom PCA
  a<-t(log(a+0.01,2))
  pc <- prcomp(a,
               center = TRUE,
               scale. = F)
  # elbow plot 
  elb<-factoextra::fviz_eig(pc)
  
  # PC1 PC2 plot
  PC12<-ggbiplot::ggbiplot(pc,
                           obs.scale = 1,
                           var.scale = 1,
                           # groups = meta[["celltype_genotype"]],
                           # colour = meta[["batch"]],
                           # fill = meta[["celltype_genotype"]],
                           # shape = meta[["batch"]],
                           ellipse = T,
                           groups = meta[["group2"]],
                           circle = F,
                           ellipse.prob = 0.68,
                           var.axes = FALSE)+
    geom_point(size = 2, aes(color=meta[["group2"]]
                             , shape=meta[["cohort"]]
    ))+
    theme(aspect.ratio=1)+
    labs(title = suffix)
  PC12<-PC12 + guides(color=guide_legend("group2")
                      ,shape=guide_legend("cohort")
  )
  
  
  pc.sample.df<-as.data.frame(pc$x) #sample
  pc.sample.df$ID<-row.names(pc.sample.df)
  
  pc.mi.df<-as.data.frame(pc$rotation)  #miRNA/mRNA
  
  # miRNA-PC1 vs PC2 colored by sample
  # PC12.sample<-ggbiplot::ggbiplot(pc,
  #                       obs.scale = 1,
  #                       var.scale = 1,
  #                       groups = meta[["celltype_genotype"]],
  #                       colour = meta[["batch"]],
  #                       ellipse = F,
  #                       circle = F,
  #                       ellipse.prob = 0.68,
  #                       var.axes = FALSE)
  
  
  print(elb)
  print(PC12)
  # print(PC12.sample)
  
  aa<-list("pc" = pc,
           "elb" = elb,
           "PC12" = PC12,
           "pc.sample.df"=pc.sample.df,
           "pc.mi.df"=pc.mi.df)
  
  return(aa)
}


# adjusted.cpm.all<-raw.to.cpm.filter(adjusted.raw.count.all, metadata.all)
# colnames(adjusted.cpm.all$df)<-metadata.all$simple_ID

# PCA.raw<-PCA_plot(cpm.all$df,metadata.all,"all_cohort")
# PCA.raw<-PCA_plot(adjusted.cpm.all$df,metadata.all,"adjusted.all_cohort")

PCA.raw<-PCA_plot(CHW.IR.cpm,CHW.IR.meta,"CHW-IR")

# flip the x and y axis and make plot with or without label
{
  PCA.raw$pc$x[,1]<- -(PCA.raw$pc$x[,1])
  PCA.raw$pc$x[,2]<- -(PCA.raw$pc$x[,2])
  PCA.raw$pc$rotation[,1]<- -(PCA.raw$pc$rotation[,1])
  PCA.raw$pc$rotation[,2]<- -(PCA.raw$pc$rotation[,2])
  meta<-CHW.IR.meta
  suffix<-"CHW-IR"
  PC12<-ggbiplot::ggbiplot(PCA.raw$pc,
                           obs.scale = 1,
                           var.scale = 1,
                           # groups = meta[["celltype_genotype"]],
                           # colour = meta[["batch"]],
                           # fill = meta[["celltype_genotype"]],
                           # shape = meta[["batch"]],
                           ellipse = T,
                           groups = meta[["group2"]],
                           circle = F,
                           ellipse.prob = 0.68,
                           var.axes = FALSE)+
    geom_point(size = 2, aes(color=meta[["group2"]]
                             , shape=meta[["cohort"]]
    ))+
    theme(aspect.ratio=1)+
    labs(title = suffix)
  PC12<-PC12 + guides(color=guide_legend("group2")
                      ,shape=guide_legend("cohort"))
  PC12
  PC12+ 
    geom_text_repel(
      aes(label = CHW.IR.meta$simple_ID),
      family = "Poppins",
      size = 3,
      min.segment.length = 0, 
      seed = 42, 
      box.padding = 0.5,
      max.overlaps = Inf,
      arrow = arrow(length = unit(0.010, "npc")),
      nudge_x = .15,
      nudge_y = .5,
      color = "grey50"
    )
  
  # eigengene table
  
  PCA.lnc.input<-PCA.raw$pc
  
  eigengene<-as.data.frame(PCA.lnc.input$rotation)
  saveRDS(eigengene,"20250610-CHW134-eigengene.RDS")
  openxlsx::write.xlsx(eigengene,"20250610-CHW134-eigengene.xlsx", rowNames = T, colNames = T)
  
}

# add label to the PCA plot
PCA.raw$PC12+ 
  geom_text_repel(
    aes(label = CHW.IR.meta$simple_ID),
    family = "Poppins",
    size = 3,
    min.segment.length = 0, 
    seed = 42, 
    box.padding = 0.5,
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = .15,
    nudge_y = .5,
    color = "grey50"
  )


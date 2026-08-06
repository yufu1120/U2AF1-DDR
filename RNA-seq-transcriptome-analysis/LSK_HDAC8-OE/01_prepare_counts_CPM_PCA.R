analysis_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
analysis_output_dir <- file.path(analysis_dir, "output")
dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)

# package
{
  library(devtools)
  library(ggbiplot)
  library(tidyverse)
  library(edgeR)
  library(ggrepel)
  library(sva)
  library(openxlsx)
setwd(analysis_dir)
  }

# setwd to raw data
getwd()
list.files()


# load raw data:

data_file_1 <- "counts/mouse_mm10/counts_ykuo_250208_mm10.txt"
data <- read.table(data_file_1, header = T, sep = "\t")
colnames(data)
data<-data[,c(6,8:ncol(data))]
row.names(data)<-data[,1]
data<-data[,2:ncol(data)]
colnames(data) <- gsub("X","",colnames(data))



metadata<-openxlsx::read.xlsx("./metadata/metadata-MIG-HDAC8-IR.xlsx")
metadata
row.names(metadata)<-metadata$Sample_ID
colnames(data)<-rownames(metadata)

all(colnames(data) %in% rownames(metadata))
index<-match(rownames(metadata),colnames(data))
data<-data[,index]

all(colnames(data) ==rownames(metadata))

raw.count.all<-data
metadata.all<-metadata

#save the metadata and raw count
# saveRDS(raw.count.all,"./counts/20250326_raw.count.all.rds")
# write.csv(raw.count.all,"./counts/20250326_raw.count.all.csv")

colnames(metadata.all)

metadata.all$cohort<-as.factor(metadata.all$cohort)

metadata.all<- metadata.all %>% 
  unite(group2, c("phenotype", "treatment"), remove=F)

# saveRDS(metadata.all, "../metadata/20250326_metadata.all.rds")



# start from here
raw.count.all<-readRDS("./counts/20250326_raw.count.all.rds")
metadata.all<-readRDS("./metadata/20250326_metadata.all.rds")


# make raw count into cpm by edgeR ================================================
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
# 
# saveRDS(cpm.all,"20250326_cpm.count.all.rds")
# write.csv(cpm.all$df,"20250326_cpm.count.all.csv")

# start from here:

cpm.all<-readRDS("./counts/20250326_cpm.count.all.rds")
metadata.all<-readRDS("./metadata/20250326_metadata.all.rds")

# split cohort:
OE.IR.cpm<-cpm.all$df
OE.IR.meta<-metadata.all

a<-OE.IR.cpm
meta<-OE.IR.meta
suffix<-"OE_IR"


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
                           # groups = meta[["group2"]],
                           colour = meta$group2,
                           fill = meta[["group2"]],
                           # shape = meta[["batch"]],
                           ellipse = T,
                           groups = meta[["group2"]],
                           circle = F,
                           ellipse.prob = 0.68,
                           var.axes = FALSE)+
    # geom_point(size = 2, aes(color=meta[["group2"]]
    #                          , shape=meta[["group2"]]
    # ))+
    theme(aspect.ratio=1)+
    labs(title = suffix)+
    geom_text_repel(
      aes(label = OE.IR.meta$Sample_ID),
      force = 1,
      direction = "both",
      nudge_x = .5,
      # nudge_y = .5,
      # family = "Poppins",
      size = 3,
      min.segment.length = 0.5, 
      seed = 42, 
      box.padding = 1,
      max.overlaps = Inf,
      arrow = arrow(length = unit(0.010, "npc")),
      
      color = "black"
    )
  PC12
  PC12<-PC12 + guides(color=guide_legend("group2"),
                      shape=guide_legend("cohort")
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



set.seed(126)
# PCA.raw<-PCA_plot(cpm.all$df,metadata.all,"all_cohort")
PCA.raw<-PCA_plot(OE.IR.cpm,OE.IR.meta,"OE_IR")


# 
# category<-unique(metadata.all$group2)


# add label to the PCA plot
p <- ggplot(pc$x, aes(PC1, PC2, colour = meta$group2)) +
  stat_ellipse(geom = "path", level=0.8, alpha=0.5)+
  # stat_ellipse(geom = "polygon", col= "black", alpha =0.1)+
  geom_text_repel(
    label = meta$Sample_ID, 
    force = 2,
    direction = "both",
    nudge_x = .5,
    nudge_y = .5,
    # family = "Poppins",
    size = 3,
    min.segment.length = 0.5, 
    seed = 42, 
    box.padding = 1.5,
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.010, "npc"))
  )+
  
  geom_point()+
  theme(aspect.ratio=1)+
  labs(title = suffix)
p

#remove sample 61362 and 61370 
setwd(analysis_dir)
cpm.all<-readRDS("./counts/20250326_cpm.count.all.rds")
metadata.all<-readRDS("./metadata/20250326_metadata.all.rds")
# split cohort:
OE.IR.cpm<-cpm.all$DEG[,c(1:4,6:12,14:16)]
OE.IR.cpm.df<-cpm.all$df[,c(1:4,6:12,14:16)]
OE.IR.meta<-metadata.all[c(1:4,6:12,14:16),]
# PCA after filtered
PCA.raw<-PCA_plot(OE.IR.cpm.df, OE.IR.meta,"OE_IR")

# draw manually but not using function to change the color:
{

  
  # PCA ==========================================================================
  a<-OE.IR.cpm.df
  meta<-OE.IR.meta
  suffix<-"OE_IR"
    # perfrom PCA
    a<-t(log(a+0.01,2))
    pc <- prcomp(a,
                 center = TRUE,
                 scale. = F)
    # elbow plot 
    elb<-factoextra::fviz_eig(pc)
    
    # PC1 PC2 plot
    group_colors <- c("WT_NIR"= "#00BFC4"      ,
                      "OEHDAC8_NIR" = "#FFA040",
                      "WT_IR" = "#000000" ,
                      "OEHDAC8_IR" = "#B856D7")  # change as needed
    
    ## light blue (MIG-NIR): #00BFC4 
    ## dark blue (MIG-IR): #0B5394
    ## orange (HDAC8-NIR): #FFA040
    ## purple (HDAC8-IR): #B856D7
    
    PC12 <- ggbiplot::ggbiplot(pc,
                               obs.scale = 1,
                               var.scale = 1,
                               groups = meta$group2,  # only needed here
                               ellipse = TRUE,
                               ellipse.prob = 0.68,
                               circle = FALSE,
                               var.axes = FALSE) +
      geom_point(aes(color = meta$group2, fill = meta$group2), size = 2, shape = 21) +
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors) +
      theme(aspect.ratio = 1) +
      labs(title = suffix, color = "Group", fill = "Group")
    PC12
    # +
    #   geom_text_repel(
    #     aes(label = OE.IR.meta$Sample_ID),
    #     force = 1,
    #     direction = "both",
    #     nudge_x = .5,
    #     # nudge_y = .5,
    #     # family = "Poppins",
    #     size = 3,
    #     min.segment.length = 0.5, 
    #     seed = 42, 
    #     box.padding = 1,
    #     max.overlaps = Inf,
    #     arrow = arrow(length = unit(0.010, "npc")),
    #     
    #     color = "black"
    #   )
    PC12
    PC12<-PC12 + guides(color=guide_legend("group2"),
                        shape=guide_legend("cohort")
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
    

  
}

# eigengene table

PCA.lnc.input<-PCA.raw$pc
eigengene<-as.data.frame(PCA.lnc.input$rotation)
saveRDS(eigengene,"20250610-MIG-eigengene.RDS")
openxlsx::write.xlsx(eigengene,"20250610-MIG-eigengene.xlsx", rowNames = T, colNames = T)

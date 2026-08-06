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
  library(openxlsx)
  
  # setwd to raw data

  setwd(analysis_dir)
  getwd()
  list.files()
}



# load raw data:

data_file_1 <- "./gene_count.xls"
data <- read.table(data_file_1, header = T, sep = "\t")
colnames(data)
row.names(data)<-data[,1]
gene_ID_table<-data[,c(1,38)]
saveRDS(gene_ID_table,'gene_ID_table.RDS')
# gene_ID_table<-readRDS('gene_ID_table.RDS')

data<-data[,c(2:37)]


metadata<-data.frame(Sample_ID = colnames(data))
metadata$celltype<-c(rep("MV4_11",12),rep("K562", 24))
metadata$cohort<-c(rep("1",12),rep("2", 24))
metadata$phenotype<-c(rep("CD531",3), rep("HDAC8", 3),
                      rep("CD531",3), rep("HDAC8", 3),
                      rep("KD_WT",3),rep("KD_K23R",3),rep("KD_K175R",3),rep("KD_4mut",3),
                      rep("KD_WT",3),rep("KD_K23R",3),rep("KD_K175R",3),rep("KD_4mut",3))
metadata$treatment<-c(rep(c("NIR","IR"),each = 6),rep(c("NIR","IR"),each = 12))
metadata$replicate<-c(rep(c(1:3),12))
metadata$group<-c(rep(LETTERS[1:12], each = 3))
metadata$simple_ID<-paste(metadata$cohort, metadata$phenotype, metadata$treatment, metadata$replicate, sep = "-")
metadata$order<-c(1:36)


row.names(metadata)<-metadata$Sample_ID
all(colnames(data) ==rownames(metadata))

all(colnames(data) %in% rownames(metadata))
index<-match(rownames(metadata),colnames(data))
data<-data[,index]

raw.count.all<-data
metadata.all<-metadata

# save the metadata and raw count
# saveRDS(raw.count.all,"./2025-04-28_raw.count.all.rds")
# write.csv(raw.count.all,"./2025-04-28_raw.count.all.csv")

colnames(metadata.all)

metadata.all$cohort<-as.factor(metadata.all$cohort)
metadata.all<- metadata.all %>% 
  unite(group2, c("phenotype", "treatment"), remove=F)


# saveRDS(metadata.all, "./2025-04-28_metadata.all.rds")


# start from here
raw.count.all<-readRDS("./counts/2025-04-28_raw.count.all.rds")
metadata.all<-readRDS("./metadata/2025-04-28_metadata.all.rds")
openxlsx::write.xlsx(metadata.all, "2025-04-28_metadata.all.xlsx")


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
# saveRDS(cpm.all, paste0(Sys.Date(),"_cpm.count.all.rds"))
# write.csv(cpm.all$df,paste0(Sys.Date(),"_cpm.count.all.csv"))
list.files()

# start from here:
setwd(analysis_dir)
cpm.all<-readRDS("./counts/2025-04-30_cpm.count.all.rds")
metadata.all<-readRDS("./metadata/2025-04-28_metadata.all.rds")





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
                           # fill = meta[["phenotype"]],
                           # shape = meta[["treatment"]],
                           ellipse = T,
                           groups = meta[["group2"]],
                           circle = F,
                           ellipse.prob = 0.68,
                           var.axes = FALSE)+
    # geom_point(size = 2, aes(color=meta[["phenotype"]]
    #                          , shape=meta[["treatment"]]
    # ))+
    theme(aspect.ratio=1)+
    labs(title = suffix)
  # +
  #   geom_text_repel(
  #     aes(label = meta$Sample_ID),
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
  
  return(aa)
}



set.seed(126)
# split cohort:
colnames(cpm.all$df)

all.cpm<-cpm.all$df
all.meta<-metadata.all
suffix<-"all"
PCA.all<-PCA_plot(OE.IR.cpm,OE.IR.meta,"all")


MV4.cpm<-cpm.all$DEG[,1:12]
MV4.cpm.df<-cpm.all$df[,1:12]
MV4.meta<-metadata.all[1:12,]
suffix<-"MV4"
PCA.MV4<-PCA_plot(MV4.cpm.df,MV4.meta,suffix)

K562.cpm<-cpm.all$DEG[,13:36]
K562.cpm.df<-cpm.all$df[,13:36]
K562.meta<-metadata.all[13:36,]
suffix<-"K562"
PCA.K562<-PCA_plot(K562.cpm.df, K562.meta, suffix)

K562.cpm<-cpm.all$df[,c(13,14,16:36)]
K562.meta<-metadata.all[c(13,14,16:36),]
suffix<-"K562-trim"
PCA.K562<-PCA_plot(K562.cpm.df, K562.meta, suffix)

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

cpm.all<-readRDS("./counts/20250326_cpm.count.all.rds")
metadata.all<-readRDS("../metadata/20250326_metadata.all.rds")
# split cohort:
OE.IR.cpm<-cpm.all$DEG[,c(1:4,6:12,14:16)]
OE.IR.cpm.df<-cpm.all$df[,c(1:4,6:12,14:16)]
OE.IR.meta<-metadata.all[c(1:4,6:12,14:16),]
# PCA after filtered
PCA.raw<-PCA_plot(OE.IR.cpm.df, OE.IR.meta,"OE_IR")

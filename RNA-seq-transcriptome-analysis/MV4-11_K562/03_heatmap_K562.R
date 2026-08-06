analysis_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
analysis_output_dir <- file.path(analysis_dir, "output")
dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)

{
  library(circlize)
  library(ggplot2)
  library(ggbiplot)
  library(colorRamp2)
  library(ComplexHeatmap)
  library(magrittr) # needed to load the pipe function '%%'
  # library("xlsx")
  library(openxlsx)
  library(dplyr)
}

# load data
setwd(analysis_dir)
cpm.all<-readRDS("./counts/2025-04-30_cpm.count.all.rds")
metadata.all<-readRDS("./metadata/2025-04-28_metadata.all.rds")

# split cohort:
# MV4.cpm<-cpm.all$DEG[,1:12]
# MV4.cpm.df<-cpm.all$df[,1:12]
# MV4.meta<-metadata.all[1:12,]
# suffix<-"MV4"
# 
# cpm2<-MV4.cpm.df

K562.cpm<-cpm.all$DEG[,c(13,14,16:26,28:36)]
K562.cpm.df<-cpm.all$df[,c(13,14,16:26,28:36)]
K562.meta<-metadata.all[c(13,14,16:26,28:36),]
suffix<-"K562"

cpm2<-K562.cpm.df

DEG.list1<-openxlsx::read.xlsx("./output/DEG/K562/2025-04-30-DEG-IvsE.xlsx", sheet =1)
DEG.list2<-openxlsx::read.xlsx("./output/DEG/K562/2025-04-30-DEG-JvsF.xlsx", sheet =1)
DEG.list3<-openxlsx::read.xlsx("./output/DEG/K562/2025-04-30-DEG-KvsG.xlsx", sheet =1)
DEG.list4<-openxlsx::read.xlsx("./output/DEG/K562/2025-04-30-DEG-LvsH.xlsx", sheet =1)

DEG1 <- DEG.list1 %>% filter(!is.na(significant))
DEG2 <- DEG.list2 %>% filter(!is.na(significant))
DEG3 <- DEG.list3 %>% filter(!is.na(significant))
DEG4 <- DEG.list4 %>% filter(!is.na(significant))

list1<-DEG1$genes 
list2<-DEG2$genes
list3<-DEG3$genes
list4<-DEG4$genes

# overlap the gene list
# Define gene list combinations
gene_lists <- list(
  com1 = list1,
  com2 = list2,
  com3 = list3,
  com4 = list4,
  com12 = union(list1, list2),
  com34 = union(list3, list4),
  com13 = union(list1, list3),
  com24 = union(list2, list4),
  com1234 = union(c(list1, list2), c(list3, list4))
)

unique(K562.meta$phenotype)

# Define right_annotation
ha <- rowAnnotation(
  phenotype = K562.meta$phenotype,
  treatment = K562.meta$treatment,
  col = list(
    phenotype = c("KD_WT" = "#DC76FF", "KD_K23R" = "#75EBFF", "KD_K175R" = "#D1FF75", "KD_4mut" = "#FFB375"),
    treatment = c("IR" = "darkred", "NIR" = "gray")
  ),
  gp = gpar(col = "black")
)


# data = cpm.count.all$df
# genelist = gene_lists[[1]]
# column_km = 4
# genelist_name = names(gene_lists)[1]

# data = cpm2
# genelist <- gene_lists[[genelist_name]]
# column_km = 3
# genelist_name = genelist_name

# Define heatmap function and with annotation of gene name
heatmap_clusterlist <- function(data, genelist, column_km = 4, genelist_name = "default") {
  cpm3 <- data[rownames(data) %in% genelist, ]
  cpm3 <- t(scale(t(cpm3))) # z-score scaling by gene
  mat  <- t(cpm3)
  
  ht <- ComplexHeatmap::Heatmap(
    mat,
    column_title = "",
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    heatmap_legend_param = list(title = "z-score", title_position = "leftcenter-rot"),
    right_annotation = ha,
    row_km = 4,
    column_km = column_km,
    row_gap = unit(1, "mm"),
    column_gap = unit(1, "mm")
  )
  
  # Save PNG
  png(file = paste0(Sys.Date(), "-heatmap-", genelist_name, "-k", column_km, ".png"),
      res = 300, units = "in", height = 6, width = 8)
  HM <- draw(ht)
  dev.off()
  
  # Extract order and clusters
  col.list <- column_order(HM)
  row.list <- row_order(HM)
  
  merge <- unlist(col.list)
  cpm_subset <- data[rownames(data) %in% genelist, ]
  cpm_reorder <- cpm_subset[merge, ]
  
  mat <- cpm_subset
  clu_df <- lapply(names(col.list), function(i) {
    data.frame(GeneID = rownames(mat[col.list[[i]], ]),
               Cluster = i,
               stringsAsFactors = FALSE)
  }) %>% do.call(rbind, .)
  
  # Export Excel
  # clu_df <- rownames_to_column(clu_df, var = "gene_name")
  clu_df<-clu_df %>%
    dplyr::left_join(.,DEG.list1[,1:3], join_by(GeneID==genes))
  write.xlsx(clu_df,
             file = paste0(Sys.Date(), "-gene-clusters-", genelist_name, "-k", column_km, ".xlsx"),
             sheetName = "Sheet1",
             col.names = TRUE, row.names = TRUE, append = FALSE)
}

getwd()
output_dir <- "./output/heatmap/K562/heatmap"
if (!dir.exists(output_dir)) dir.create(output_dir)
setwd(output_dir)

# Loop over all gene list combinations
for (genelist_name in names(gene_lists)) {
  message("Processing: ", genelist_name)
  genelist <- gene_lists[[genelist_name]]
  heatmap_clusterlist(data = cpm2, genelist = genelist, column_km = 3, genelist_name = genelist_name)
}

for (genelist_name in names(gene_lists)) {
  message("Processing: ", genelist_name)
  genelist <- gene_lists[[genelist_name]]
  heatmap_clusterlist(data = cpm2, genelist = genelist, column_km = 4, genelist_name = genelist_name)
}


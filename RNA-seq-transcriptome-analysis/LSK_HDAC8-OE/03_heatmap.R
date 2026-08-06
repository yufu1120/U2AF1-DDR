analysis_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
analysis_output_dir <- file.path(analysis_dir, "output")
dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)

{
library(circlize)
library(ggplot2)
library(ggbiplot)
library(colorRamp2)
library(ComplexHeatmap)
library(magrittr)
library(xlsx)
}

# load data
{
  setwd(analysis_dir)
  cpm.count.all<-readRDS("./counts/20250326_cpm.count.all.rds")
  metadata.all<-readRDS("./metadata/20250326_metadata.all.rds")

  
  cpm2<-cpm.count.all$df[,c(1:4,6:12,14:16)]
  metadata.all<-metadata.all[c(1:4,6:12,14:16),]
  
  
  DEG1<-read.table("./output/sig.BvsA.csv", sep = ",", header = T)
  DEG2<-read.table("./output/sig.DvsC.csv", sep = ",", header = T)
  DEG3<-read.table("./output/sig.CvsA.csv", sep = ",", header = T)
  DEG4<-read.table("./output/sig.DvsB.csv", sep = ",", header = T)
  
  list1<-DEG1$genes 
  list2<-DEG2$genes
  list3<-DEG3$genes
  list4<-DEG4$genes
  
}

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


# Define right_annotation
ha <- rowAnnotation(
  phenotype = metadata.all$phenotype,
  treatment = metadata.all$treatment,
  col = list(
    phenotype = c("OEHDAC8" = "#F8766D", "WT" = "#00BFC4"),
    treatment = c("IR" = "darkred", "NIR" = "gray")
  ),
  gp = gpar(col = "black")
)


# data = cpm.count.all$df
# genelist = gene_lists[[1]]
# column_km = 4
# genelist_name = names(gene_lists)[1]



# Define heatmap function
heatmap_clusterlist <- function(data, genelist, column_km = 4, genelist_name = "default") {
  cpm3 <- data[rownames(data) %in% genelist, ]
  cpm3 <- t(scale(t(cpm3))) # z-score scaling by gene
  mat <- t(cpm3)
  
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
  write.xlsx(clu_df,
             file = paste0(Sys.Date(), "-gene-clusters-", genelist_name, "-k", column_km, ".xlsx"),
             sheetName = "Sheet1",
             col.names = TRUE, row.names = TRUE, append = FALSE)
}

getwd()
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

# draw a single heatmap
gene_lists<- readxl::read_xlsx("./2024-05-28-gene-clusters5-com34-k3.xlsx")
gene_list<-gene_lists[gene_lists$Cluster==3,]
# message("Processing: ", genelist_name)
# genelist <- gene_lists[[genelist_name]]

genelist<-gene_list$GeneID
genelist_name<-gene_list
# heatmap_clusterlist(data = cpm2, genelist = genelist, column_km = 4, genelist_name = genelist_name)

data <- cpm2
genelist <- genelist
column_km = 3
genelist_name = "CM_to_MIG"
  
  {
  cpm3 <- data[rownames(data) %in% genelist, ]
  cpm3 <- t(scale(t(cpm3))) # z-score scaling by gene
  mat <- t(cpm3)
  
  ht <- ComplexHeatmap::Heatmap(
    mat,
    column_title = "",
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    cluster_rows = F,
    # cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    heatmap_legend_param = list(title = "z-score", title_position = "leftcenter-rot"),
    right_annotation = ha,
    row_km = 0,
    column_km = column_km,
    row_gap = unit(1, "mm"),
    column_gap = unit(1, "mm")
  )
  ht
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
  write.xlsx(clu_df,
             file = paste0(Sys.Date(), "-gene-clusters-", genelist_name, "-k", column_km, ".xlsx"),
             sheetName = "Sheet1",
             col.names = TRUE, row.names = TRUE, append = FALSE)
}

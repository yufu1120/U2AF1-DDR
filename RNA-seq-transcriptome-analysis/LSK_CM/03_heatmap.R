analysis_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
analysis_output_dir <- file.path(analysis_dir, "output")
dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)

{
library(circlize)
library(ggplot2)
library(ggbiplot)
# install.packages("colorRamp2")
library(colorRamp2)
library(ComplexHeatmap)
}

# load data
setwd(analysis_dir)
cpm.all<-readRDS("./counts/20230425_cpm.count.all.rds")
metadata.all<-readRDS("./metadata/20230425_metadata.all.rds")

adjusted.cpm.all<-cpm.all

# split cohort:
cm_count_columns <- c(1, 3:9, 11:16)
CHW.IR.cpm<-cpm.all$df[,cm_count_columns]
CHW.IR.meta<-metadata.all
stopifnot(ncol(CHW.IR.cpm) == nrow(CHW.IR.meta))

adjusted.cpm.all$df<-CHW.IR.cpm


# draw heatmap for top 100 variavle genes
{
sd <- apply(adjusted.cpm.all$df,1,sd)
means <- rowMeans(adjusted.cpm.all$df)
cv<-sd/means
hist(cv)

cv.ordered <- order(cv,decreasing=T)
var_top_100 <- adjusted.cpm.all$df[cv.ordered[1:500],]


# ComplexHeatmap::Heatmap(adjusted.cpm.all)


col_fun = colorRamp2(c(-1, 0, 1), c("darkblue", "white", "yellow"))
# col_fun(seq(-3, 3))



var_top_100<-log(var_top_100+min_offset,2)
var_top_100<-t(scale(t(var_top_100))) #performed by column

ComplexHeatmap::Heatmap(var_top_100,
                        column_title  = "Top 100 most variable genes",
                        col = col_fun,
                        row_names_gp = gpar(fontsize = 0),
                        column_names_gp  = gpar(fontsize = 10),
                        cluster_columns = T,
                        show_row_names = F,
                        heatmap_legend_param = list(title = "row-normalized z-score",
                                                    title_position = "leftcenter-rot")
)
}

# draw heatmap for WT, CHW, NIR, IR:
setwd(analysis_output_dir)
list.files()

# load data
CHW_WT_DEG<-read.table("sig.BvsA.csv", sep = ",", header = T)
CHW_WT_IR_DEG<-read.table("sig.DvsC.csv", sep = ",", header = T)
WT_IR_DEG<-read.table("sig.CvsA.csv", sep = ",", header = T)
CHW_IR_DEG<-read.table("sig.DvsB.csv", sep = ",", header = T)

list1<-CHW_WT_DEG$genes
list2<-CHW_WT_IR_DEG$genes
list3<-WT_IR_DEG$genes
list4<-CHW_IR_DEG$genes

# union_list<-list1
# union_list<-list3
union_list<-union(list1, list2)
union_list34<-union(list3, list4)
union_list<-union(list1, list3)
union_list2<-union(union_list,list1)
union_list3<-union(union_list2,list2)

union_list<-union_list2
union_list<-union_list3
# cpm
getwd()
cpm.count.all<-readRDS(file.path(analysis_dir, "counts", "20230425_cpm.count.all.rds"))

cpm.count.all$df<-cpm.count.all$df[,cm_count_columns]

# # 1. get top CV genes
# sd <- apply(cpm.count.all$df,1,sd)
# means <- rowMeans(cpm.count.all$df)
# cv<-sd/means
# hist(cv)
# cv.ordered <- order(cv,decreasing=T)
# var_top_100 <- cpm.count.all$df[cv.ordered[1:500],]
# cpm2<-log(var_top_100+0.01,2)


# 2.get gene list from DEGs
min_offset<-min(cpm.count.all$df[cpm.count.all$df!=0])
cpm2<-log(cpm.count.all$df+min_offset,2)

# cpm2<-cpm2[row.names(cpm2) %in% union_list,]
colnames(cpm2)
dim(cpm2) # 2354 gene x 12 samples
cpm2<-t(scale(t(cpm2))) #scaling performed by column (gene)


heatmap_filter(cpm2, union_list)



# col_fun = colorRamp2(c(-1, 0, 1), c("darkblue", "white", "yellow"))
head(cpm.count.all$df)
mat<-cpm2
ComplexHeatmap::Heatmap(mat,
                        column_title  = "",
                        # col = col_fun,
                        row_names_gp = gpar(fontsize = 0),
                        column_names_gp  = gpar(fontsize = 10),
                        cluster_columns = T,
                        show_row_names = F,
                        heatmap_legend_param = list(title = "z-score",
                                                    title_position = "leftcenter-rot")
)


# to make heatmap by DEGs from each group and combine...

list1<-CHW_WT_DEG$genes
list2<-CHW_WT_IR_DEG$genes
list3<-WT_IR_DEG$genes
list4<-CHW_IR_DEG$genes

union_list34<-union(list3, list4)

# define group annotation
{
CHW.IR.meta$phenotype[CHW.IR.meta$phenotype=="WT"]<-'Ctrl'
CHW.IR.meta$phenotype[CHW.IR.meta$phenotype=="CHW"]<-'CM'
col_letters = c("Ctrl" = "gray", "CM" = "red")
ha = rowAnnotation(phenotype = CHW.IR.meta$phenotype,
                   treatment = CHW.IR.meta$treatment,
                   col = list(phenotype = c("CM" = "#F8766D", "Ctrl" = "#00BFC4"),
                              treatment = c("IR" = "darkred", "NIR" = "gray")),
                   gp = gpar(col = "black"))
draw(ha)
}

heatmap_filter<-function(data, genelist, column_km){
  cpm3<-cpm2[row.names(cpm2) %in% genelist,]
  cpm3<-t(scale(t(cpm3))) #scaling performed by column (gene)
  
  mat<-t(cpm3)
  ComplexHeatmap::Heatmap(mat,
                          column_title  = "",
                          # col = col_fun,
                          row_names_gp = gpar(fontsize = 10),
                          column_names_gp  = gpar(fontsize = 0),
                          cluster_rows = T,
                          cluster_columns = T,
                          show_row_names = T,
                          show_column_names = F,
                          heatmap_legend_param = list(title = "z-score",
                                                      title_position = "leftcenter-rot"),
                          right_annotation = ha,
                          row_km = 4,
                          column_km =column_km,
                          row_gap = unit(1, "mm"),
                          column_gap = unit(1, "mm"),

  )
}
# 
# ht1<-heatmap_filter(cpm2, list1)
# ht2<-heatmap_filter(cpm2, list2)
# ht3<-heatmap_filter(cpm2, list3)
# ht4<-heatmap_filter(cpm2, list4)



png(file=paste0(Sys.Date(),"-heatmap-com34-k3.png"), res=300, units="in", height=6, width=8)
set.seed(123)
ht34<-heatmap_filter(cpm2, union_list34, column_km=3)
HM <-draw(ht34)
graphics.off()


# extract the order from cluster 
col.dend <- column_dend(HM)  #If needed, extract row dendrogram
col.list <- column_order(HM)  #Extract clusters (output is a list)
row.list <- row_order(HM)

lapply(col.list, function(x) length(x))  #check/confirm size gene clusters
merge<-c(col.list$`3`, col.list$`2`, col.list$`1`)

cpm34<-cpm2[row.names(cpm2) %in% union_list34,]
cpm34_reorder<-cpm34[merge,]


library(magrittr) # needed to load the pipe function '%%'
library("xlsx")
mat<-cpm34
clu_df <- lapply(names(col.list), function(i){
  out <- data.frame(GeneID = rownames(mat[col.list[[i]],]),
                    Cluster = i,
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)

table(clu_df$Cluster)

#export
write.xlsx(clu_df, file = paste0(Sys.Date(),"-gene-clusters5-com34-k3.xlsx"),
           sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)



######## draw heatmap for CM vs Ctrl (NIR) and cm vs Ctrl (IR) #########
setwd(analysis_output_dir)
list.files()

# load data
CHW_WT_DEG<-read.table("sig.BvsA.csv", sep = ",", header = T)
CHW_WT_IR_DEG<-read.table("sig.DvsC.csv", sep = ",", header = T)
# WT_IR_DEG<-read.table("sig.CvsA.csv", sep = ",", header = T)
# CHW_IR_DEG<-read.table("sig.DvsB.csv", sep = ",", header = T)

list1<-CHW_WT_DEG$genes
list2<-CHW_WT_IR_DEG$genes
# list3<-WT_IR_DEG$genes
# list4<-CHW_IR_DEG$genes

# union_list<-list1
# union_list<-list3
union_list12<-union(list1, list2)
# union_list34<-union(list3, list4)
# union_list<-union(list1, list3)

union_list<-union_list12


# cpm
getwd()
cpm.count.all<-readRDS(file.path(analysis_dir, "counts", "20230425_cpm.count.all.rds"))

cpm.count.all$df<-cpm.count.all$df[,cm_count_columns]

# # 1. get top CV genes
# sd <- apply(cpm.count.all$df,1,sd)
# means <- rowMeans(cpm.count.all$df)
# cv<-sd/means
# hist(cv)
# cv.ordered <- order(cv,decreasing=T)
# var_top_100 <- cpm.count.all$df[cv.ordered[1:500],]
# cpm2<-log(var_top_100+0.01,2)


# 2.get gene list from DEGs
min_offset<-min(cpm.count.all$df[cpm.count.all$df!=0])
cpm2<-log(cpm.count.all$df+min_offset,2)

# cpm2<-cpm2[row.names(cpm2) %in% union_list,]
colnames(cpm2)
dim(cpm2) # 2354 gene x 12 samples
cpm2<-t(scale(t(cpm2))) #scaling performed by column (gene)





# col_fun = colorRamp2(c(-1, 0, 1), c("darkblue", "white", "yellow"))
head(cpm.count.all$df)
mat<-cpm2
ComplexHeatmap::Heatmap(mat,
                        column_title  = "",
                        # col = col_fun,
                        row_names_gp = gpar(fontsize = 0),
                        column_names_gp  = gpar(fontsize = 10),
                        cluster_columns = T,
                        show_row_names = F,
                        heatmap_legend_param = list(title = "z-score",
                                                    title_position = "leftcenter-rot")
)

# define group annotation
{
  CHW.IR.meta$phenotype[CHW.IR.meta$phenotype=="WT"]<-'Ctrl'
  CHW.IR.meta$phenotype[CHW.IR.meta$phenotype=="CHW"]<-'CM'
  col_letters = c("Ctrl" = "gray", "CM" = "red")
  ha = rowAnnotation(phenotype = CHW.IR.meta$phenotype,
                     treatment = CHW.IR.meta$treatment,
                     col = list(phenotype = c("CM" = "#F8766D", "Ctrl" = "#00BFC4"),
                                treatment = c("IR" = "darkred", "NIR" = "gray")),
                     gp = gpar(col = "black"))
  draw(ha)
}

heatmap_filter<-function(data, genelist, column_km){
  cpm3<-cpm2[row.names(cpm2) %in% genelist,]
  cpm3<-t(scale(t(cpm3))) #scaling performed by column (gene)
  
  mat<-t(cpm3)
  ComplexHeatmap::Heatmap(mat,
                          column_title  = "",
                          # col = col_fun,
                          row_names_gp = gpar(fontsize = 10),
                          column_names_gp  = gpar(fontsize = 0),
                          cluster_rows = T,
                          cluster_columns = T,
                          show_row_names = T,
                          show_column_names = F,
                          heatmap_legend_param = list(title = "z-score",
                                                      title_position = "leftcenter-rot"),
                          right_annotation = ha,
                          row_km = 4,
                          column_km =column_km,
                          row_gap = unit(1, "mm"),
                          column_gap = unit(1, "mm"),
                          
  )
}
# 
# ht1<-heatmap_filter(cpm2, list1)
# ht2<-heatmap_filter(cpm2, list2)
# ht3<-heatmap_filter(cpm2, list3)
# ht4<-heatmap_filter(cpm2, list4)

union_list12<-union(list1, list2)

png(file=paste0(Sys.Date(),"-heatmap-com12-k2.png"), res=300, units="in", height=6, width=8)
set.seed(123)
ht12<-heatmap_filter(cpm2, union_list12, column_km=2)
HM <-draw(ht12)
graphics.off()


# extract the order from cluster 
col.dend <- column_dend(HM)  #If needed, extract row dendrogram
col.list <- column_order(HM)  #Extract clusters (output is a list)
row.list <- row_order(HM)

lapply(col.list, function(x) length(x))  #check/confirm size gene clusters
merge<-c(col.list$`3`, col.list$`2`, col.list$`1`)

cpm12<-cpm2[row.names(cpm2) %in% union_list12,]
cpm12_reorder<-cpm12[merge,]


library(magrittr) # needed to load the pipe function '%%'
library("xlsx")
mat<-cpm12
clu_df <- lapply(names(col.list), function(i){
  out <- data.frame(GeneID = rownames(mat[col.list[[i]],]),
                    Cluster = i,
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)

table(clu_df$Cluster)

#export
write.xlsx(clu_df, file = paste0(Sys.Date(),"-gene-clusters2-com12.xlsx"),
           sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

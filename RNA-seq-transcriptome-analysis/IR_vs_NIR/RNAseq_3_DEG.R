

# package
{
library(dplyr)
library(tidyr)
library(edgeR)
library(cmapR)
library(ggplot2)

library(clusterProfiler)
library(KEGGREST)
library(gprofiler2)

library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(wesanderson)
library(viridis)
library(EnhancedVolcano)
library(data.table)
library(scales)
library(gprofiler2)
library(ReactomePA)
library(stringr)
library(xlsx)
library(ggh4x)
}


# setwd
setwd(r"(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\HDAC8)")
# setwd("~/HDAC8/20230906-CHW-134")

# load data
cpm.all<-readRDS("./counts/20230425_cpm.count.all.rds")
metadata.all<-readRDS("./metadata/20230425_metadata.all.rds")

# 20230906 only have CHW-1,3,4 and WT-1,3,4
CHW.IR.cpm<-cpm.all$DEG[,c(1,3:5,7:9,11:13,15:16)]
CHW.IR.cpm.df<-cpm.all$df[,c(1,3:5,7:9,11:13,15:16)]
CHW.IR.meta<-metadata.all[c(1,3:5,7:9,11:13,15:16),]


######################################### cohort: CHW ==========================
{
y<-CHW.IR.cpm
metadata<-CHW.IR.meta
metadata$group<-as.factor(metadata$group)

colnames(metadata)

# design matrix
row.names(y$samples)<-metadata$simple_ID
# gender<-metadata$gender
group<-metadata$group
# group2<-metadata$group2

# design <- model.matrix(~0+group+gender)  # blocking gender
design <- model.matrix(~0+group)
rownames(design) <- row.names(y$samples)

y <- estimateDisp(y,design,robust=TRUE)
fit <- glmQLFit(y,design,robust=TRUE)

# start to perform indicated comparison:


setwd(paste0(r"(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\HDAC8)","\\20230906-CHW-134\\output"))

date<-Sys.Date()
cohort<-"CHW"
e1<-"groupA"
e2<-"groupB"
e3<-"groupC"
e4<-"groupD"

con <- makeContrasts(com1=eval(parse(text =e2))-eval(parse(text =e1)),
                     com2=eval(parse(text =e4))-eval(parse(text =e3)),
                     com3=eval(parse(text =e3))-eval(parse(text =e1)),
                     com4=eval(parse(text =e4))-eval(parse(text =e2)),
                     levels = design)


for (i in 1:length(colnames(con))){
  assign(paste0(cohort,".qlf.",colnames(con)[i]), glmQLFTest(fit, contrast=con[,i]))
  assign(paste0(cohort,".qlf.toptags.",colnames(con)[i]),
         as.data.frame(topTags(eval(parse(text = paste0(cohort,".qlf.",colnames(con)[i]))), n=nrow(y))))
  a<-eval(parse(text=paste0(cohort,".qlf.toptags.",colnames(con)[i])))
  colnames(a)
  a<-a %>%
    arrange(desc(logFC))
  assign(paste0(cohort,".qlf.toptags.",colnames(con)[i]),a)
  
  # write.csv(eval(parse(text=paste0(cohort,".qlf.toptags.",colnames(con)[i]))),
  #                  paste(date,cohort,colnames(con)[i],"DEG.csv",sep = "-"))
  }
}

       
# =filter out the DEGs =============================================================================
# testing chunk
group1<-"A"
group2<-"B"
qlf.toptags<-CHW.qlf.toptags.com1
cpmtable<-CHW.IR.cpm.df
metatable<-CHW.IR.meta

# this function will keep only the DEGs in the table
DEGsfilter<-function(group1,group2,qlf.toptags){
df<-tibble::rownames_to_column(cpm.all$df)
sample.index<-which(metadata.all$group %in% c(group1,group2))
mat<-df[,c(1, sample.index+1)]

df<-left_join(qlf.toptags, mat, by=c("genes" = "rowname"))

df1 <- df %>%
  # 1st filtering: filter out those low-expression genes () 
  # filter(rowSums(df[,c(7:ncol(df))]) > 1) %>%
  # 2nd filtering: p-val <= 0.05
  filter( PValue <= 0.05) %>%
  # 3rd filtering: logFC >= log2(1) or logFC <= log2(1/2)
  filter( logFC >= log2(1.5) | logFC <= log2(0.67) ) %>%
  # sort by fold change
  arrange(desc(logFC),desc(PValue))
write.csv(df1,paste0("sig.",group2,"vs",group1,".csv"))

num.up<-sum(df1[["logFC"]]>= log2(1.5))
num.down<-sum(df1[["logFC"]]<= log2(0.67))
num.sum=num.up+num.down
message(
  paste0("There are ",num.up," up-regulated genes and ",num.down," down-regulated genes, in total ", num.sum, " genes, in comparison of group ",group2," to ",group1,"."))

return(df1)

}


sig.BvsA<-DEGsfilter("A","B",CHW.qlf.toptags.com1)
sig.DvsC<-DEGsfilter("C","D",CHW.qlf.toptags.com2)
sig.CvsA<-DEGsfilter("A","C",CHW.qlf.toptags.com3)
sig.DvsB<-DEGsfilter("B","D",CHW.qlf.toptags.com4)

# this function will keep all the genes but labeling the DEGs
DEGsfilter2<-function(group1,group2,qlf.toptags,cpmtable, metatable){
  date<-Sys.Date()
  df<-tibble::rownames_to_column(cpmtable)
  sample.index<-which(metatable$group %in% c(group1,group2))
  mat<-df[,c(1, sample.index+1)]
  
  df<-left_join(qlf.toptags, mat, by=c("genes" = "rowname"))
  
  index.sig.up<-df$PValue <= 0.05 & df$logFC >= log2(2)
  df$significant[index.sig.up]<-"up"
  
  index.sig.down<-df$PValue <= 0.05 & df$logFC <= log2(0.5)
  df$significant[index.sig.down]<-"down"
  
  write.xlsx(df,paste0(date,"-DEG-",group2,"vs",group1,".xlsx"))
  
  num.up<-sum(index.sig.up)
  num.down<-sum(index.sig.down)
  num.sum=num.up+num.down
  message(
    paste0("There are ",num.up," up-regulated genes and ",num.down," down-regulated genes, in total ", num.sum, " genes, in comparison of group ",group2," to ",group1,"."))
  
  return(df)
  
}
options(java.parameters = "-Xmx8000m")
sig.BvsA<-DEGsfilter2("A","B",CHW.qlf.toptags.com1, CHW.IR.cpm.df, CHW.IR.meta)
sig.DvsC<-DEGsfilter2("C","D",CHW.qlf.toptags.com2, CHW.IR.cpm.df, CHW.IR.meta)
sig.CvsA<-DEGsfilter2("A","C",CHW.qlf.toptags.com3, CHW.IR.cpm.df, CHW.IR.meta)
sig.DvsB<-DEGsfilter2("B","D",CHW.qlf.toptags.com4, CHW.IR.cpm.df, CHW.IR.meta)


# this function will keep all the genes but labeling the DEGs and using threshold of 1.5/0.67
DEGsfilter_v2_1.5<-function(group1,group2,qlf.toptags,cpmtable, metatable){
  date<-Sys.Date()
  df<-tibble::rownames_to_column(cpmtable)
  sample.index<-which(metatable$group %in% c(group1,group2))
  mat<-df[,c(1, sample.index+1)]
  
  df<-left_join(qlf.toptags, mat, by=c("genes" = "rowname"))
  
  index.sig.up<-df$PValue <= 0.05 & df$logFC >= log2(1.5)
  df$significant[index.sig.up]<-"up"
  
  index.sig.down<-df$PValue <= 0.05 & df$logFC <= log2(0.67)
  df$significant[index.sig.down]<-"down"
  
  write.xlsx(df,paste0(date,"-DEG-",group2,"vs",group1,".xlsx"))
  
  num.up<-sum(index.sig.up)
  num.down<-sum(index.sig.down)
  num.sum=num.up+num.down
  message(
    paste0("There are ",num.up," up-regulated genes and ",num.down," down-regulated genes, in total ", num.sum, " genes, in comparison of group ",group2," to ",group1,"."))
  
  return(df)
  
}
options(java.parameters = "-Xmx8000m")
sig.BvsA<-DEGsfilter_v2_1.5("A","B",CHW.qlf.toptags.com1, CHW.IR.cpm.df, CHW.IR.meta)
sig.DvsC<-DEGsfilter_v2_1.5("C","D",CHW.qlf.toptags.com2, CHW.IR.cpm.df, CHW.IR.meta)
sig.CvsA<-DEGsfilter_v2_1.5("A","C",CHW.qlf.toptags.com3, CHW.IR.cpm.df, CHW.IR.meta)
sig.DvsB<-DEGsfilter_v2_1.5("B","D",CHW.qlf.toptags.com4, CHW.IR.cpm.df, CHW.IR.meta)

# =volcano plot =============================================================================
# Use EnhancedVolcano package
# result<-qlf.toptags.BvsA

draw.colcano.plot<-function(result, title, col){
  FCcutoff<-1
  dum<-0.01
  pCutoff<-0.05
  p<-EnhancedVolcano(result, 
                  # lab=result$genes,  #decide to have label or not
                  lab=NA,
                  # selectLab = result$genes[(result$logFC<= -2|result$logFC>= 2 )&result$PValue <= dum],
                  x="logFC", 
                  y="PValue", 
                  FCcutoff = FCcutoff, 
                  pCutoff = pCutoff,
                  xlim=c(-7,7),
                  ylim=c(0,5),
                  pointSize = 1.2, 
                  labSize = 3, 
                  colAlpha = 0.8, 
                  title=title, 
                  subtitle=NULL, 
                  caption=NULL,  
                  drawConnectors = TRUE,
                  # maxoverlapsConnectors = 20,
                  # legendPosition="bottom",
                  legendPosition = "none"
                  # col= col #Ns,Log2FC.p-value, FC and p-value
                  # legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                  # 'p-value & Log (base 2) FC'),
                  # legendLabSize = 10
  ) 
  
  png(file=paste0(Sys.Date(),"-volcano-",title,"-.png"), res=300, units="in", height=3.6, width=6)
  print(p)
  graphics.off()
  return(p)
}
library(patchwork)
title<-"CHW-NIR vs WT-NIR"
p1<-draw.colcano.plot(CHW.qlf.toptags.com1,"CHW-NIR vs WT-NIR")
p2<-draw.colcano.plot(CHW.qlf.toptags.com2,"CHW-IR vs WT-IR")
p3<-draw.colcano.plot(CHW.qlf.toptags.com3,"Ctrl-IR vs Ctrl-NIR", col=c("gray", "gray", "gray", "#00BFC4"))
p4<-draw.colcano.plot(CHW.qlf.toptags.com4,"CM-IR vs CM-NIR", col=c("gray", "gray", "gray", "#F8766D")) #Ns,Log2FC.p-value, FC and p-value)

getwd()
p1<-draw.colcano.plot(CHW.qlf.toptags.com1,"CHW-NIR vs WT-NIR")
p2<-draw.colcano.plot(CHW.qlf.toptags.com2,"CHW-IR vs WT-IR")
p3<-draw.colcano.plot(CHW.qlf.toptags.com3,"WT-IR vs WT-NIR")
p4<-draw.colcano.plot(sig.DvsB,"CHW-IR vs CHW-NIR")


sig.BvsA<-DEGsfilter("A","B",CHW.qlf.toptags.com1)
sig.DvsC<-DEGsfilter("C","D",CHW.qlf.toptags.com2)
sig.CvsA<-DEGsfilter("A","C",CHW.qlf.toptags.com3)
sig.DvsB<-DEGsfilter("B","D",CHW.qlf.toptags.com4)

(p1+p2)/(p3+p4)
# (p5+p6)/(p7+p8)
(p9+p10)/(p11+p12)

# =GSEA =============================================================================
# using the expression "fold change" * "p-value" as matrix 

# df<-CHW.qlf.com3
# comparison="WT(IR vs NIR)"

# set the comparison group
cohort<-c("CHW.qlf.com1", "CHW.qlf.com2", "CHW.qlf.com3", "CHW.qlf.com4")
cohort_tag<-c("CHW-NIR vs WT-NIR", "CHW-IR vs WT-IR", "WT-IR vs WT-NIR", "CHW-IR vs CHW-NIR")

# hallmark
category<-"H"
m_db = msigdbr(species = "Mus musculus", category = category) #"Mus musculus" or "Homo sapiens"
m_db_2 <- m_db[,c("gs_name", "gene_symbol")]

# load the hallmark gene set:
category<-"H"
file<-r"(C:\Users\yufu\OneDrive - City of Hope National Medical Center\rotation\1_Dr. Kuo\U2AF1\RNAseq\GSEA(X)\geneset\mice\mh.all.v0.3.symbols.gmt)"
no_col <- max(count.fields(file, sep = "\t"))
hallmark <- read.table(file,sep="\t", fill=TRUE,col.names=c("gene","link",1:(no_col-2)))

hallmark<-hallmark[,-2]
colnames(hallmark)[1]<-"gs_name"
hallmark<-hallmark %>%
  pivot_longer(
    cols = !gs_name,
    names_to = "X",
    values_to = "gene_symbol"
  ) %>%
  drop_na() %>%
  .[.$gene_symbol!="",-2]
m_db_3 <-hallmark 

# m_db_3%>%
#   group_by(gs_name) %>%
#   count()

setwd(r'(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\HDAC8\20230906-CHW-134)')
getwd()
for ( k in 1: length(cohort)){
  
  df<-get(cohort[k])
  comparison<-cohort_tag[k]
  
# calculate the matrix
{
A<-df$table$logFC
B<--log(df$table$PValue,10)
score <-A*B
# score <-A
names(score) <- row.names(df$table)
score <- sort(score, decreasing = T)
# score
}

length(unique(row.names(df$table)))
length(row.names(df$table))

# write.table(score, paste0(date,"_GSEA_",category,"_",comparison,"-score.txt"), sep = "\t")

# run GSEA

em <- GSEA(score, TERM2GENE=m_db_3, pvalueCutoff=1, pAdjustMethod="BH", eps=0)
write.csv(em, paste0(".\\output\\GSEA\\GSEA-output-table\\",Sys.Date(),"_GSEA_",category,"_",comparison,".csv"))

em.sig<-em@result[em@result$p.adjust<=1,]
# em.sig.h<-em.sig[em.sig$p.adjust<=0.05,]

# enrichment plots

# em.sig.plot<-em.sig$ID[em.sig$p.adjust<=0.05]
em.sig.plot<-em.sig$ID

for (i in 1:length(em.sig.plot)){
  geneSetID<-em.sig.plot[i]
  NES<-em@result$NES[em@result$ID==geneSetID]
  NES<-format(NES, digits=4)
  pvalue<-em@result$p.adjust[em@result$ID==geneSetID]
  pvalue<-format(pvalue, digits=4)
  FDR<-em@result$qvalues[em@result$ID==geneSetID]
  FDR<-format(FDR, digits=4)
  p1<-gseaplot2(em, geneSetID=geneSetID, title=geneSetID)+
    theme(plot.title = element_text(size=0.5))
    # annotate("text", x=0.95, y=0.9, label= paste("NES=",NES), hjust = 1)+
    # annotate("text", x=0.95, y=0.85, label= paste("adj-p-value=",pvalue), hjust = 1)+
    # annotate("text", x=0.95, y=0.8, label= paste("FDR=",FDR), hjust = 1)
  p1
  ggsave(paste0(".\\output\\GSEA\\GSEA-plot\\",Sys.Date(),"_GSEA_",category,"_",comparison,"-",str_sub(geneSetID, 10, str_count(geneSetID)),".png"),dpi = 300)
  
}

# bobble plot 1 (depricated):
# em.sig %>%
#   slice_max(n=50, order_by = NES) %>%
#   arrange(NES)%>%
#   mutate(abs.NES=abs(NES))%>%
#   mutate(Description=factor(Description, levels=Description))%>%
#   ggplot(aes(x=NES, y=Description, colour=p.adjust, size=qvalues))+
#   geom_point()+
#   ggtitle(paste0("GSEA_",category,"_",comparison))+
#   # xlim(-5, -12)+
#   labs(x="NES", y="Hallmark", colour="adjusted p-value", size="qvalues")+
#   scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")),breaks = c(0.05,0.2,0.4,0.6,0.8))

}

# bobble plot 2 : plot for hallmark top 15 and last 15
setwd(r'(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\HDAC8\20230906-CHW-134)')
plots<-".\\output\\GSEA\\GSEA-output-table\\"
group<-c("WT-IR vs WT-NIR","CHW-NIR vs WT-NIR","CHW-IR vs WT-IR","CHW-IR vs CHW-NIR")
getwd()
k=1
for (i in 1:length(group)){
  em<-read.csv(paste0(".\\output\\GSEA\\GSEA-output-table\\2024-04-12_GSEA_H_",group[i],".csv"), header = T, sep = ",")
  
  
  em.plot<-em %>%
    mutate(abs.NES=abs(NES))%>%
    dplyr::arrange(NES) %>%
    slice(c(1:15, 36:50)) %>%
    mutate("ID_short"=unlist(lapply(ID, function(x) paste(strsplit(x,"_")[[1]][-1],collapse=" ")))) %>%
    
    mutate(ID_short=factor(ID_short, levels=ID_short))
  
  png(paste(plots,group[i],"-hmk-top15.png", sep=""), height=6, width=8, res=300, units="in")
  
  p<-em.plot %>%
    # arrange(desc(NES))%>%
    ggplot(aes(x=NES, y=ID_short, size=abs.NES, colour=pvalue))+
    geom_point()+
    # ggtitle(paste0("GSEA_",category,"_",comparison))+
    xlim(-2, 2)+
    labs(x="NES", y="", colour="p-value", size="absolute NES")+
    # scale_size_continuous(limits=c(0,1))+
    # scale_color_gradientn(colours = rainbow(5))+
    scale_color_gradientn(colours = rev(wes_palette(n=5, name="Zissou1")),breaks = c(0.05,0.2,0.4,0.6,0.8))+
    # scale_color_brewer(palette = "Dark2")+
    # scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")),breaks = c(0.05,0.2,0.4,0.6,0.8))+
    theme_bw()+
    theme(panel.background = element_rect(fill = "white"),
          # axis.text = element_text(size = 6, face="bold", colour = "black"),
          axis.text = element_text(color="black"),
    )+
    # guides(colour = guide_legend(order = 1),size = guide_legend(order = 2))+
    
    labs(title = group[i])
  
  print(p)
  graphics.off()
}


# C2
category<-"C2"
m_db = msigdbr(species = "Mus musculus", category = category) #"Mus musculus" or "Homo sapiens"
m_db_2 <- m_db[,c("gs_name", "gene_symbol")]
em <- GSEA(score, TERM2GENE=m_db_2, pvalueCutoff=1, pAdjustMethod="BH", eps=0)
write.csv(em, paste0(date,"_GSEA_",category,"_",comparison,".csv"))


# em.sig<-em@result[em@result$p.adjust<=1,]
em.sig<-em@result[em@result$p.adjust<=0.05,]
em.sig.C2<-em.sig[em.sig$p.adjust<=0.05,]
# write.csv(em.sig, "20230615_GSEA_C2_CDDvsWT.sig.csv")
colnames(em.sig)
em.sig %>%
  slice_max(n=50, order_by = NES) %>%
  arrange(NES)%>%
  mutate(abs.NES=abs(NES))%>%
  mutate(Description=factor(Description, levels=Description))%>%
  ggplot(aes(x=NES, y=Description, colour=p.adjust, size=qvalues))+
  geom_point()+
  ggtitle(paste0("GSEA_",category,"_",comparison))+
  xlim(-2, 2)+
  labs(x="NES", y="C2", colour="adjusted p-value", size="qvalues")+
  scale_colour_stepsn(colours=c("#00BFC4",wes_palette(name="BottleRocket2")[1:7]),
                      breaks = c(0,0.05,0.2,0.5,1),
                      limits = c(0,1),
                      guide = guide_coloursteps(even.steps = T,
                                                show.limits = T))
em.sig %>%
  slice_max(n=50, order_by = NES) %>%
  arrange(NES)%>%
  mutate(abs.NES=abs(NES))%>%
  mutate(Description=factor(Description, levels=Description))%>%
  ggplot(aes(x=NES, y=Description, colour=p.adjust, size=qvalues))+
  geom_point()+
  ggtitle(paste0("GSEA_",category,"_",comparison))+
  xlim(-2, 2)+
  labs(x="NES", y="C2", colour="adjusted p-value", size="qvalues")+
  scale_colour_stepsn(colours=c("#F8766D",wes_palette(name="BottleRocket2")[1:7]),
                      breaks = c(0,0.05,0.2,0.5,1),
                      limits = c(0,1),
                      guide = guide_coloursteps(even.steps = T,
                                                show.limits = T))



em.sig %>%
  slice_max(n=50, order_by = -NES) %>%
  arrange(NES)%>%
  mutate(abs.NES=abs(NES))%>%
  mutate(Description=factor(Description, levels=Description))%>%
  ggplot(aes(x=NES, y=Description, colour=p.adjust, size=qvalues))+
  geom_point()+
  ggtitle(paste0("GSEA_",category,"_",comparison))+
  # xlim(-5, -12)+
  labs(x="NES", y="C2", colour="adjusted p-value", size="qvalues")+
  scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")))

em.sig.plot<-em.sig$ID[em.sig$p.adjust<=0.05]
for (i in 1:length(em.sig.plot)){
  geneSetID<-em.sig.plot[i]
  NES<-em@result$NES[em@result$ID==geneSetID]
  NES<-format(NES, digits=4)
  pvalue<-em@result$p.adjust[em@result$ID==geneSetID]
  pvalue<-format(pvalue, digits=4)
  FDR<-em@result$qvalues[em@result$ID==geneSetID]
  FDR<-format(FDR, digits=4)
  p1<-gseaplot2(em, geneSetID=geneSetID, title=geneSetID)+
    theme(plot.title = element_text(size=0.5))+
    annotate("text", x=0.95, y=0.9, label= paste("NES=",NES), hjust = 1)+
    annotate("text", x=0.95, y=0.85, label= paste("adj-p-value=",pvalue), hjust = 1)+
    annotate("text", x=0.95, y=0.8, label= paste("FDR=",FDR), hjust = 1)
  p1
  ggsave(paste0(date,"_GSEA_",category,"_",comparison,"-",str_sub(geneSetID, 10, str_count(geneSetID)),".png"),dpi = 300)
  
}

# C5
category<-"C5"
subcategory<-"BP"
m_db = msigdbr(species = "Mus musculus", category = category, subcategory=subcategory) #"Mus musculus" or "Homo sapiens"
m_db_2 <- m_db[,c("gs_name", "gene_symbol")]
em <- GSEA(score, TERM2GENE=m_db_2, pvalueCutoff=1, pAdjustMethod="BH", eps=0)
write.csv(em, paste0(date,"_GSEA_",category,"_",comparison,".csv"))

# em.sig<-em@result[em@result$p.adjust<=1,]
em.sig<-em@result[em@result$p.adjust<=0.5,]
em.sig.C5<-em.sig[em.sig$p.adjust<=0.05,]
colnames(em.sig)
em.sig %>%
  slice_max(n=50, order_by = NES) %>%
  arrange(NES)%>%
  mutate(abs.NES=abs(NES))%>%
  mutate(Description=factor(Description, levels=Description))%>%
  ggplot(aes(x=NES, y=Description, colour=p.adjust, size=abs.NES))+
  geom_point()+
  ggtitle(paste0("GSEA_",category,"_",subcategory,"_",comparison))+
  # xlim(-5, -12)+
  labs(x="NES", y="C5", colour="adjusted p-value", size="qvalues")+
  scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")))

em.sig %>%
  slice_max(n=50, order_by = -NES) %>%
  arrange(NES)%>%
  mutate(abs.NES=abs(NES))%>%
  mutate(Description=factor(Description, levels=Description))%>%
  ggplot(aes(x=NES, y=Description, colour=p.adjust, size=abs.NES))+
  geom_point()+
  ggtitle(paste0("GSEA_",category,"_",subcategory,"_",comparison))+
  # xlim(-5, -12)+
  labs(x="NES", y="C5", colour="adjusted p-value", size="qvalues")+
  scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")))

em.sig.plot<-em.sig$ID[em.sig$p.adjust<=0.05]

for (i in 1:length(em.sig.plot)){
  geneSetID<-em.sig.plot[i]
  NES<-em@result$NES[em@result$ID==geneSetID]
  NES<-format(NES, digits=4)
  pvalue<-em@result$p.adjust[em@result$ID==geneSetID]
  pvalue<-format(pvalue, digits=4)
  FDR<-em@result$qvalues[em@result$ID==geneSetID]
  FDR<-format(FDR, digits=4)
  p1<-gseaplot2(em, geneSetID=geneSetID, title=geneSetID)+
    theme(plot.title = element_text(size=0.5))+
    annotate("text", x=0.95, y=0.9, label= paste("NES=",NES), hjust = 1)+
    annotate("text", x=0.95, y=0.85, label= paste("adj-p-value=",pvalue), hjust = 1)+
    annotate("text", x=0.95, y=0.8, label= paste("FDR=",FDR), hjust = 1)
  p1
  ggsave(paste0(date,"_GSEA_",category,"_",comparison,"-",str_sub(geneSetID, 10, str_count(geneSetID)),".png"),dpi = 300)
  
}


# BRCA1

# load the BRCA1 gene set:
# category<-"BRCA1"
# file<-r"(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\HDAC8\GSEA\geneset\BRCA1.txt)"
# no_col <- max(count.fields(file, sep = "\t"))
# BRCA <- read.table(file,sep="\t",fill=TRUE,col.names=c("gene","link",1:(no_col-2)))
# 
# # BRCA<-read.table(file, header = F,fill = TRUE, )
# BRCA<-BRCA[,-2]
# 
# colnames(BRCA)[1]<-"gs_name"
# BRCA<-BRCA %>% 
#   pivot_longer(
#     cols = !gs_name, 
#     names_to = "X", 
#     values_to = "gene_symbol"
#   ) %>%
#   drop_na() %>%
#   .[.$gene_symbol!="",-2]
# m_db_3 <-BRCA


em <- GSEA(score, TERM2GENE=m_db_3, pvalueCutoff=1, pAdjustMethod="BH", eps=0)
write.csv(em, paste0(date,"_GSEA_",category,"_",comparison,".csv"))

# em.sig<-em@result[em@result$p.adjust<=1,]
em.sig<-em@result[em@result$p.adjust<=0.5,]

# em.sig.WT<-em.sig
# em.sig.CHW<-em.sig

# em.sig.BRCA1<-em.sig[em.sig$p.adjust<=0.05,]
# colnames(em.sig.WT)


em.sig<-em.sig.WT
em.sig<-em.sig.CHW

n_up<-sum(em.sig$NES>0)
n_down<-sum(em.sig$NES<0)

em.sig %>%
  slice_max(n=n_up, order_by = NES) %>%
  arrange(NES)%>%
  mutate(abs.NES=abs(NES))%>%
  mutate(Description=factor(Description, levels=Description))%>%
  ggplot(aes(x=NES, y=Description, colour=p.adjust, size=abs.NES))+
  geom_point()+
  ggtitle(paste0("GSEA_",category,"_",subcategory,"_",comparison))+
  # xlim(-5, -12)+
  labs(x="NES", y="pathway", colour="adjusted p-value", size="abs.NES")+
  # scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")))
  scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")),breaks = c(0.05,0.2,0.4,0.6,0.8))
em.sig %>%
  slice_max(n=n_down, order_by = -NES) %>%
  arrange(-NES)%>%
  mutate(abs.NES=abs(NES))%>%
  mutate(Description=factor(Description, levels=Description))%>%
  ggplot(aes(x=NES, y=Description, colour=p.adjust, size=abs.NES))+
  geom_point()+
  ggtitle(paste0("GSEA_",category,"_",subcategory,"_",comparison))+
  # xlim(-5, -12)+
  labs(x="NES", y="pathway", colour="adjusted p-value", size="abs.NES")+
  scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")))
  # scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")),breaks = c(0.05,0.2,0.4,0.6,0.8))



em.sig.plot<-


em.sig.plot<-em.sig$ID[em.sig$p.adjust<=0.05]
for (i in 1:length(em.sig.plot)){
  geneSetID<-em.sig.plot[i]
  NES<-em@result$NES[em@result$ID==geneSetID]
  NES<-format(NES, digits=4)
  pvalue<-em@result$p.adjust[em@result$ID==geneSetID]
  pvalue<-format(pvalue, digits=4)
  FDR<-em@result$qvalues[em@result$ID==geneSetID]
  FDR<-format(FDR, digits=4)
  p1<-gseaplot2(em, geneSetID=geneSetID, title=geneSetID)+
    theme(plot.title = element_text(size=0.5))+
    annotate("text", x=0.95, y=0.9, label= paste("NES=",NES), hjust = 1)+
    annotate("text", x=0.95, y=0.85, label= paste("adj-p-value=",pvalue), hjust = 1)+
    annotate("text", x=0.95, y=0.8, label= paste("FDR=",FDR), hjust = 1)
  p1
  ggsave(paste0(date,"_GSEA_",category,"_",comparison,"-",str_sub(geneSetID, 10, str_count(geneSetID)),".png"),dpi = 300)
  
}


# BRCA1

# load the BRCA1 gene set:

# download the gmt file from GSEA website and change the extension to .txt
# category<-"BRCA1"
# file<-r"(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\HDAC8\GSEA\geneset\BRCA1.txt)"

category<-"BRCA1-v2"
file<-r"(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\HDAC8\GSEA\geneset\BRCA.v2023.2.Mm.txt)"


no_col <- max(count.fields(file, sep = "\t"))
BRCA <- read.table(file,sep="\t",fill=TRUE,col.names=c("gene","link",1:(no_col-2)))

# BRCA<-read.table(file, header = F,fill = TRUE, )
BRCA<-BRCA[,-2]

colnames(BRCA)[1]<-"gs_name"
BRCA<-BRCA %>%
  pivot_longer(
    cols = !gs_name,
    names_to = "X",
    values_to = "gene_symbol"
  ) %>%
  drop_na() %>%
  .[.$gene_symbol!="",-2]
m_db_3 <-BRCA

em <- GSEA(score, TERM2GENE=m_db_3, pvalueCutoff=1, pAdjustMethod="BH", eps=0, minGSSize = 5)
write.csv(em, paste0(date,"_GSEA_",category,"_",comparison,".csv"))



#Create GSEA Plot for significant pathways
# devtools::install_github("nicolash2/gggsea")

# =try fgsea =============================================================================

# 
library(fgsea)
library(data.table)
library(ggplot2)

data(examplePathways)
data(exampleRanks)
set.seed(42)
em <- GSEA(score, TERM2GENE=m_db_2, pvalueCutoff=1, pAdjustMethod="BH", eps=0)

# Load the pathways into a named list
pathways.hallmark <- fgsea::gmtPathways(r"(C:\Users\yufu\OneDrive - City of Hope National Medical Center\rotation\1_Dr. Kuo\U2AF1\RNAseq\GSEA\geneset\mice\mh.all.v0.3.symbols.gmt)")

# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes <- fgsea(pathways = pathways.hallmark, 
                  stats    = score,
                  minSize  = 15,
                  maxSize  = 500)
head(fgseaRes[order(pval), ])


# short conclusion:
# clusterprofiler and the fgsea will have the same results (and both uses the prerank list),
# but different from java version(even using the prerank list) 
# java used weighted_p2 will be little bit closer
analysis_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
analysis_output_dir <- file.path(analysis_dir, "output")
dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)


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


setwd(analysis_dir)

# load data
cpm.all<-readRDS("./counts/20250326_cpm.count.all.rds")
metadata.all<-readRDS("./metadata/20250326_metadata.all.rds")

OE.IR.cpm<-cpm.all$DEG[,c(1:4,6:12,14:16)]
OE.IR.cpm.df<-cpm.all$df[,c(1:4,6:12,14:16)]
OE.IR.meta<-metadata.all[c(1:4,6:12,14:16),]
}

######################################### cohort: HDAC8-OE ==========================
{
y<-OE.IR.cpm
metadata<-OE.IR.meta
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
if (!dir.exists("./output")) {
  dir.create("./output")
}
setwd("./output")

date<-Sys.Date()
cohort<-"HDAC8_OE"
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


sig.BvsA<-DEGsfilter("A","B",HDAC8_OE.qlf.toptags.com1)
sig.DvsC<-DEGsfilter("C","D",HDAC8_OE.qlf.toptags.com2)
sig.CvsA<-DEGsfilter("A","C",HDAC8_OE.qlf.toptags.com3)
sig.DvsB<-DEGsfilter("B","D",HDAC8_OE.qlf.toptags.com4)

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
sig.BvsA<-DEGsfilter_v2_1.5("A","B",HDAC8_OE.qlf.toptags.com1, OE.IR.cpm.df, OE.IR.meta)
  sig.DvsC<-DEGsfilter_v2_1.5("C","D",HDAC8_OE.qlf.toptags.com2, OE.IR.cpm.df, OE.IR.meta)
sig.CvsA<-DEGsfilter_v2_1.5("A","C",HDAC8_OE.qlf.toptags.com3, OE.IR.cpm.df, OE.IR.meta)
sig.DvsB<-DEGsfilter_v2_1.5("B","D",HDAC8_OE.qlf.toptags.com4, OE.IR.cpm.df, OE.IR.meta)


# =GSEA =============================================================================
# 
# using the expression "fold change" * "p-value" as matrix 
{
  setwd(analysis_output_dir)
  if (!dir.exists("./GSEA")) {
    dir.create("./GSEA")
  }
  setwd("./GSEA")
  
  
  # set the comparison group
  cohort<-c("HDAC8_OE.qlf.com1", "HDAC8_OE.qlf.com2", "HDAC8_OE.qlf.com3", "HDAC8_OE.qlf.com4")
  cohort_tag<-c("OE-NIR vs MIG-NIR", "OE-IR vs MIG-IR", "MIG-IR vs MIG-NIR", "OE-IR vs OE-NIR")
  
  # test trunk
  # k=4
  # df<-get(cohort[k])
  # comparison<-cohort_tag[k]
  k=1
  for ( k in 1: length(cohort)){
    setwd(file.path(analysis_output_dir, "GSEA"))
    # df<-get(cohort[k])
    comparison<-cohort_tag[k]
    saved_date<-"2025-08-06"
    if (!dir.exists(paste0("./",k,"-",comparison))) {
      dir.create(paste0("./",k,"-",comparison))
    }
    setwd(paste0("./",k,"-",comparison))
    
    # calculate the matrix
    # {
    #   A<-df$table$logFC
    #   B<--log(df$table$PValue,10)
    #   score <-A*B
    #   # score <-A
    #   names(score) <- row.names(df$table)
    #   score <- sort(score, decreasing = T)
    #   # score
    # }
    # write.table(score, paste0(date,"_GSEA_",comparison,"-score.txt"), sep = "\t")
    # saveRDS(score, paste0(date,"_GSEA_",comparison,"-score.RDS"))
    # or read the score RDS:
    # score <- readRDS(paste0(saved_date,"_GSEA_",comparison,"-score.RDS"))
        # run GSEA
    {

      # (1) hallmark
      # {
      #   category<-"H"
      #   m_db = msigdbr(species = "Mus musculus", category = category) #"Mus musculus" or "Homo sapiens"
      #   m_db_2 <- m_db[,c("gs_name", "gene_symbol")]
      # 
      #   # # instaed, load the hallmark gene set from local:
      #   # category<-"H"
      #   # file<-"./mh.all.v0.3.symbols.gmt"
      #   # no_col <- max(count.fields(file, sep = "\t"))
      #   # hallmark <- read.table(file,sep="\t", fill=TRUE,col.names=c("gene","link",1:(no_col-2)))
      #   #
      #   # hallmark<-hallmark[,-2]
      #   # colnames(hallmark)[1]<-"gs_name"
      #   # hallmark<-hallmark %>%
      #   #   pivot_longer(
      #   #     cols = !gs_name,
      #   #     names_to = "X",
      #   #     values_to = "gene_symbol"
      #   #   ) %>%
      #   #   drop_na() %>%
      #   #   .[.$gene_symbol!="",-2]
      #   # m_db_3 <-hallmark
      #   #
      # 
      #   em <- GSEA(score, TERM2GENE=m_db_2, pvalueCutoff=1, pAdjustMethod="BH", eps=0)
      #   write.csv(em, paste0("./",Sys.Date(),"_GSEA_",category,"_",comparison,".csv"))
      #   saveRDS(em, paste0("./",Sys.Date(),"_GSEA_",category,"_",comparison,".RDS"))
      #   em<-readRDS(paste0("./",saved_date,"_GSEA_",category,"_",comparison,".RDS"))
      #   em.sig<-em@result[em@result$p.adjust<=1,]
      #   
      #   # plotting
      #   # {
      #   #   # enrichment plots
      #   #   {
      #   #     em.sig.plot<-em.sig$ID
      #   #     for (i in 1:length(em.sig.plot)){
      #   #       geneSetID<-em.sig.plot[i]
      #   #       NES<-em@result$NES[em@result$ID==geneSetID]
      #   #       NES<-format(NES, digits=4)
      #   #       pvalue<-em@result$p.adjust[em@result$ID==geneSetID]
      #   #       pvalue<-format(pvalue, digits=4)
      #   #       FDR<-em@result$qvalues[em@result$ID==geneSetID]
      #   #       FDR<-format(FDR, digits=4)
      #   #       p1<-gseaplot2(em, geneSetID=geneSetID, title=geneSetID)+
      #   #         theme(plot.title = element_text(size=0.5))
      #   #       # annotate("text", x=0.95, y=0.9, label= paste("NES=",NES), hjust = 1)+
      #   #       # annotate("text", x=0.95, y=0.85, label= paste("adj-p-value=",pvalue), hjust = 1)+
      #   #       # annotate("text", x=0.95, y=0.8, label= paste("FDR=",FDR), hjust = 1)
      #   #       p1
      #   #       ggsave(paste0("./",Sys.Date(),"_GSEA_",category,"_",comparison,"-",str_sub(geneSetID, 10, str_count(geneSetID)),".png"),dpi = 300)
      #   #     }
      #   #   }
      #   #   
      #   #   
      #   #   # bobble plot
      #   #   {
      #   #     # bobble plot 1 (depricated):
      #   #     {
      #   #       # em.sig %>%
      #   #       #   slice_max(n=50, order_by = NES) %>%
      #   #       #   arrange(NES)%>%
      #   #       #   mutate(abs.NES=abs(NES))%>%
      #   #       #   mutate(Description=factor(Description, levels=Description))%>%
      #   #       #   ggplot(aes(x=NES, y=Description, colour=p.adjust, size=qvalues))+
      #   #       #   geom_point()+
      #   #       #   ggtitle(paste0("GSEA_",category,"_",comparison))+
      #   #       #   # xlim(-5, -12)+
      #   #       #   labs(x="NES", y="Hallmark", colour="adjusted p-value", size="qvalues")+
      #   #       #   scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")),breaks = c(0.05,0.2,0.4,0.6,0.8))
      #   #     }
      #   #     
      #   #     # bobble plot 2 : plot for hallmark top 15 and last 15
      #   #     {
      #   #       em<-read.csv(paste0("./", Sys.Date() ,"_GSEA_H_",cohort_tag[k],".csv"), header = T, sep = ",")
      #   #       em.plot<-em %>%
      #   #         mutate(abs.NES=abs(NES))%>%
      #   #         dplyr::arrange(NES) %>%
      #   #         slice(c(1:15, 36:50)) %>%
      #   #         mutate("ID_short"=unlist(lapply(ID, function(x) paste(strsplit(x,"_")[[1]][-1],collapse=" ")))) %>%
      #   #         mutate(ID_short=factor(ID_short, levels=ID_short))
      #   #       
      #   #       png(paste("./",Sys.Date(),"-bobble-",cohort_tag[k],"-hmk-top15.png", sep=""), height=6, width=8, res=300, units="in")
      #   #       p<-em.plot %>%
      #   #         # arrange(desc(NES))%>%
      #   #         ggplot(aes(x=NES, y=ID_short, size=abs.NES, colour=pvalue))+
      #   #         geom_point()+
      #   #         # ggtitle(paste0("GSEA_",category,"_",comparison))+
      #   #         xlim(-2.5, 2.5)+
      #   #         labs(x="NES", y="", colour="p-value", size="absolute NES")+
      #   #         # scale_size_continuous(limits=c(0,1))+
      #   #         # scale_color_gradientn(colours = rainbow(5))+
      #   #         scale_color_gradientn(colours = rev(wes_palette(n=5, name="Zissou1")),breaks = c(0.05,0.2,0.4,0.6,0.8))+
      #   #         # scale_color_brewer(palette = "Dark2")+
      #   #         # scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")),breaks = c(0.05,0.2,0.4,0.6,0.8))+
      #   #         theme_bw()+
      #   #         theme(panel.background = element_rect(fill = "white"),
      #   #               # axis.text = element_text(size = 6, face="bold", colour = "black"),
      #   #               axis.text = element_text(color="black"),
      #   #         )+
      #   #         # guides(colour = guide_legend(order = 1),size = guide_legend(order = 2))+
      #   #         
      #   #         labs(title = cohort_tag[k])
      #   #       
      #   #       print(p)
      #   #       graphics.off()
      #   #     }
      #   #   }
      #   # }
      # 
      #   
      # 
      #   
      # 
      # }




      # C2
      # {
      #   category<-"C2"
      #   m_db = msigdbr(species = "Mus musculus", category = category) #"Mus musculus" or "Homo sapiens"
      #   m_db_2 <- m_db[,c("gs_name", "gene_symbol")]
      #   
      #   
      #   em <- GSEA(score, TERM2GENE=m_db_2, pvalueCutoff=1, pAdjustMethod="BH", eps=0)
      #   # write.csv(em, paste0("./",Sys.Date(),"_GSEA_",category,"_",comparison,".csv"))
      #   # saveRDS(em, paste0("./",Sys.Date(),"_GSEA_",category,"_",comparison,".RDS"))
      #   em<-readRDS(paste0("./",saved_date,"_GSEA_",category,"_",comparison,".RDS"))
      #   em.sig<-em@result[em@result$p.adjust<=1,]
      #   # em.sig<-em@result[em@result$p.adjust<=0.05,]
      #   
      #   # plotting
      #   # {
      #   #   
      #   #   
      #   #   {
      #   #     png(paste("./",Sys.Date(),"-bobble-",cohort_tag[k],"-",category,"-top25.png", sep=""), height=8, width=10, res=300, units="in")
      #   #     p<-em.sig %>%
      #   #       arrange(desc(NES)) %>%  # Sort NES in descending order
      #   #       slice_head(n = 25) %>%  # Select top 25
      #   #       bind_rows(
      #   #         em.sig %>%
      #   #           arrange(NES) %>%  # Sort NES in ascending order
      #   #           slice_head(n = 25)  # Select bottom 25
      #   #       ) %>%
      #   #       arrange(NES) %>%  # Ensure correct order for plotting
      #   #       mutate(abs.NES = abs(NES)) %>%
      #   #       mutate(Description = factor(Description, levels = Description)) %>%
      #   #       ggplot(aes(x = NES, y = Description, colour = p.adjust, size = qvalues)) +
      #   #       geom_point() +
      #   #       ggtitle(paste0("GSEA_", category, "_", comparison)) +
      #   #       xlim(-2.5, 2.5) +
      #   #       labs(x = "NES", y = "C2", colour = "adjusted p-value", size = "qvalues") +
      #   #       scale_colour_stepsn(colours = c("#F8766D", wes_palette(name = "BottleRocket2")[1:7]),
      #   #                           breaks = c(0, 0.05, 0.2, 0.5, 1),
      #   #                           limits = c(0, 1),
      #   #                           guide = guide_coloursteps(even.steps = TRUE,
      #   #                                                     show.limits = TRUE))
      #   #     print(p)
      #   #     graphics.off()
      #   #   }
      #   #   
      #   #   
      #   #   # em.sig %>%
      #   #   #   slice_max(n=50, order_by = NES) %>%
      #   #   #   arrange(NES)%>%
      #   #   #   mutate(abs.NES=abs(NES))%>%
      #   #   #   mutate(Description=factor(Description, levels=Description))%>%
      #   #   #   ggplot(aes(x=NES, y=Description, colour=p.adjust, size=qvalues))+
      #   #   #   geom_point()+
      #   #   #   ggtitle(paste0("GSEA_",category,"_",comparison))+
      #   #   #   xlim(-2, 2)+
      #   #   #   labs(x="NES", y="C2", colour="adjusted p-value", size="qvalues")+
      #   #   #   scale_colour_stepsn(colours=c("#00BFC4",wes_palette(name="BottleRocket2")[1:7]),
      #   #   #                       breaks = c(0,0.05,0.2,0.5,1),
      #   #   #                       limits = c(0,1),
      #   #   #                       guide = guide_coloursteps(even.steps = T,
      #   #   #                                                 show.limits = T))
      #   #   # em.sig %>%
      #   #   #   slice_max(n=50, order_by = NES) %>%
      #   #   #   arrange(NES)%>%
      #   #   #   mutate(abs.NES=abs(NES))%>%
      #   #   #   mutate(Description=factor(Description, levels=Description))%>%
      #   #   #   ggplot(aes(x=NES, y=Description, colour=p.adjust, size=qvalues))+
      #   #   #   geom_point()+
      #   #   #   ggtitle(paste0("GSEA_",category,"_",comparison))+
      #   #   #   xlim(-2, 2)+
      #   #   #   labs(x="NES", y="C2", colour="adjusted p-value", size="qvalues")+
      #   #   #   scale_colour_stepsn(colours=c("#F8766D",wes_palette(name="BottleRocket2")[1:7]),
      #   #   #                       breaks = c(0,0.05,0.2,0.5,1),
      #   #   #                       limits = c(0,1),
      #   #   #                       guide = guide_coloursteps(even.steps = T,
      #   #   #                                                 show.limits = T))
      #   #   #
      #   #   #
      #   #   #
      #   #   # em.sig %>%
      #   #   #   slice_max(n=50, order_by = -NES) %>%
      #   #   #   arrange(NES)%>%
      #   #   #   mutate(abs.NES=abs(NES))%>%
      #   #   #   mutate(Description=factor(Description, levels=Description))%>%
      #   #   #   ggplot(aes(x=NES, y=Description, colour=p.adjust, size=qvalues))+
      #   #   #   geom_point()+
      #   #   #   ggtitle(paste0("GSEA_",category,"_",comparison))+
      #   #   #   # xlim(-5, -12)+
      #   #   #   labs(x="NES", y="C2", colour="adjusted p-value", size="qvalues")+
      #   #   #   scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")))
      #   #   
      #   #   em.sig.plot<-em.sig$ID[em.sig$p.adjust<=0.05]
      #   #   for (i in 1:length(em.sig.plot)){
      #   #     geneSetID<-em.sig.plot[i]
      #   #     NES<-em@result$NES[em@result$ID==geneSetID]
      #   #     NES<-format(NES, digits=4)
      #   #     pvalue<-em@result$p.adjust[em@result$ID==geneSetID]
      #   #     pvalue<-format(pvalue, digits=4)
      #   #     FDR<-em@result$qvalues[em@result$ID==geneSetID]
      #   #     FDR<-format(FDR, digits=4)
      #   #     p1<-gseaplot2(em, geneSetID=geneSetID, title=geneSetID)+
      #   #       theme(plot.title = element_text(size=0.5))+
      #   #       annotate("text", x=0.95, y=0.9, label= paste("NES=",NES), hjust = 1)+
      #   #       annotate("text", x=0.95, y=0.85, label= paste("adj-p-value=",pvalue), hjust = 1)+
      #   #       annotate("text", x=0.95, y=0.8, label= paste("FDR=",FDR), hjust = 1)
      #   #     p1
      #   #     ggsave(paste0(date,"_GSEA_",category,"_",comparison,"-",str_sub(geneSetID, 10, str_count(geneSetID)),".png"),dpi = 300)
      #   #   }
      #   # }
      #   
      # }


      # C5
      {
        category<-"C5"
        subcategory<-"BP"
        # m_db = msigdbr(species = "Mus musculus", category = category, subcategory=subcategory) #"Mus musculus" or "Homo sapiens"
        # m_db_2 <- m_db[,c("gs_name", "gene_symbol")]
        # em <- GSEA(score, TERM2GENE=m_db_2, pvalueCutoff=1, pAdjustMethod="BH", eps=0)
        # write.csv(em, paste0("./",Sys.Date(),"_GSEA_",category,"_",comparison,".csv"))
        # saveRDS(em, paste0("./",Sys.Date(),"_GSEA_",category,"_",comparison,".RDS"))
        em<-readRDS(paste0("./",saved_date,"_GSEA_",category,"_",comparison,".RDS"))
        # em.sig<-em@result[em@result$p.adjust<=1,]
        
        # plotting
        {

          # {
          #   png(paste("./",Sys.Date(),"-bobble-",cohort_tag[k],"-",category,"-top25.png", sep=""), height=8, width=10, res=300, units="in")
          # 
          #   p<-em.sig %>%
          #     arrange(desc(NES)) %>%  # Sort NES in descending order
          #     slice_head(n = 25) %>%  # Select top 25
          #     bind_rows(
          #       em.sig %>%
          #         arrange(NES) %>%  # Sort NES in ascending order
          #         slice_head(n = 25)  # Select bottom 25
          #     ) %>%
          #     arrange(NES) %>%  # Ensure correct order for plotting
          #     mutate(abs.NES = abs(NES)) %>%
          #     mutate(Description = factor(Description, levels = Description)) %>%
          #     ggplot(aes(x = NES, y = Description, colour = p.adjust, size = qvalues)) +
          #     geom_point() +
          #     ggtitle(paste0("GSEA_", category, "_", comparison)) +
          #     xlim(-2.5, 2.5) +
          #     labs(x = "NES", y = "C5", colour = "adjusted p-value", size = "qvalues") +
          #     scale_colour_stepsn(colours = c("#F8766D", wes_palette(name = "BottleRocket2")[1:7]),
          #                         breaks = c(0, 0.05, 0.2, 0.5, 1),
          #                         limits = c(0, 1),
          #                         guide = guide_coloursteps(even.steps = TRUE,
          #                                                   show.limits = TRUE))
          #   print(p)
          #   graphics.off()
          # }

          # em.sig.plot<-em.sig$ID[em.sig$p.adjust<=0.05]
        
          gene_list<-openxlsx::read.xlsx("../GSEA-summary-bobble_backup.xlsx")
          em.sig.plot<-gene_list$ID
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
            ggsave(paste0(Sys.Date(),"_GSEA_",category,"_",comparison,"-",geneSetID,".png"),
                   dpi = 300,
                   width = 1940, height = 1434,
                   units = "px")
          }
        }

      }



      # BRCA1
      # {
      #   # load the BRCA1 gene set:
      #   # category<-"BRCA1"
      #   # file<-"./BRCA1.txt"
      #   # no_col <- max(count.fields(file, sep = "\t"))
      #   # BRCA <- read.table(file,sep="\t",fill=TRUE,col.names=c("gene","link",1:(no_col-2)))
      #   #
      #   # # BRCA<-read.table(file, header = F,fill = TRUE, )
      #   # BRCA<-BRCA[,-2]
      #   #
      #   # colnames(BRCA)[1]<-"gs_name"
      #   # BRCA<-BRCA %>%
      #   #   pivot_longer(
      #   #     cols = !gs_name,
      #   #     names_to = "X",
      #   #     values_to = "gene_symbol"
      #   #   ) %>%
      #   #   drop_na() %>%
      #   #   .[.$gene_symbol!="",-2]
      #   # m_db_3 <-BRCA
      #   
      #   
      #   em <- GSEA(score, TERM2GENE=m_db_3, pvalueCutoff=1, pAdjustMethod="BH", eps=0)
      #   write.csv(em, paste0("./",Sys.Date(),"_GSEA_",category,"_",comparison,".csv"))
      #   saveRDS(em, paste0("./",Sys.Date(),"_GSEA_",category,"_",comparison,".RDS"))
      #   # em.sig<-em@result[em@result$p.adjust<=1,]
      #   em.sig<-em@result[em@result$p.adjust<=0.5,]
      #   
      #   # em.sig.WT<-em.sig
      #   # em.sig.CHW<-em.sig
      #   
      #   # em.sig.BRCA1<-em.sig[em.sig$p.adjust<=0.05,]
      #   # colnames(em.sig.WT)
      #   
      #   
      #   em.sig<-em.sig.WT
      #   em.sig<-em.sig.CHW
      #   
      #   n_up<-sum(em.sig$NES>0)
      #   n_down<-sum(em.sig$NES<0)
      #   
      #   em.sig %>%
      #     slice_max(n=n_up, order_by = NES) %>%
      #     arrange(NES)%>%
      #     mutate(abs.NES=abs(NES))%>%
      #     mutate(Description=factor(Description, levels=Description))%>%
      #     ggplot(aes(x=NES, y=Description, colour=p.adjust, size=abs.NES))+
      #     geom_point()+
      #     ggtitle(paste0("GSEA_",category,"_",subcategory,"_",comparison))+
      #     # xlim(-5, -12)+
      #     labs(x="NES", y="pathway", colour="adjusted p-value", size="abs.NES")+
      #     # scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")))
      #     scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")),breaks = c(0.05,0.2,0.4,0.6,0.8))
      #   em.sig %>%
      #     slice_max(n=n_down, order_by = -NES) %>%
      #     arrange(-NES)%>%
      #     mutate(abs.NES=abs(NES))%>%
      #     mutate(Description=factor(Description, levels=Description))%>%
      #     ggplot(aes(x=NES, y=Description, colour=p.adjust, size=abs.NES))+
      #     geom_point()+
      #     ggtitle(paste0("GSEA_",category,"_",subcategory,"_",comparison))+
      #     # xlim(-5, -12)+
      #     labs(x="NES", y="pathway", colour="adjusted p-value", size="abs.NES")+
      #     scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")))
      #   # scale_colour_stepsn(colours=rev(wes_palette(n=5, name="Zissou1")),breaks = c(0.05,0.2,0.4,0.6,0.8))
      #   
      #   
      #   
      #   
      #   
      #   em.sig.plot<-em.sig$ID[em.sig$p.adjust<=0.05]
      #   for (i in 1:length(em.sig.plot)){
      #     geneSetID<-em.sig.plot[i]
      #     NES<-em@result$NES[em@result$ID==geneSetID]
      #     NES<-format(NES, digits=4)
      #     pvalue<-em@result$p.adjust[em@result$ID==geneSetID]
      #     pvalue<-format(pvalue, digits=4)
      #     FDR<-em@result$qvalues[em@result$ID==geneSetID]
      #     FDR<-format(FDR, digits=4)
      #     p1<-gseaplot2(em, geneSetID=geneSetID, title=geneSetID)+
      #       theme(plot.title = element_text(size=0.5))+
      #       annotate("text", x=0.95, y=0.9, label= paste("NES=",NES), hjust = 1)+
      #       annotate("text", x=0.95, y=0.85, label= paste("adj-p-value=",pvalue), hjust = 1)+
      #       annotate("text", x=0.95, y=0.8, label= paste("FDR=",FDR), hjust = 1)
      #     p1
      #     ggsave(paste0(date,"_GSEA_",category,"_",comparison,"-",str_sub(geneSetID, 10, str_count(geneSetID)),".png"),dpi = 300)
      #     
      #   }
      # }

      


      # BRCA1-v2
      # {
      # 
      #   # load the BRCA1 gene set:
      #   # download the gmt file from GSEA website and change the extension to .txt
      #   # category<-"BRCA1"
      #   # file<-"./BRCA1.txt"
      # 
      #   category<-"BRCA1-v2"
      #   file<-"./BRCA.v2023.2.Mm.txt"
      #   no_col <- max(count.fields(file, sep = "\t"))
      #   BRCA <- read.table(file,sep="\t",fill=TRUE,col.names=c("gene","link",1:(no_col-2)))
      #   BRCA<-BRCA[,-2]
      #   colnames(BRCA)[1]<-"gs_name"
      #   BRCA<-BRCA %>%
      #     pivot_longer(
      #       cols = !gs_name,
      #       names_to = "X",
      #       values_to = "gene_symbol"
      #     ) %>%
      #     drop_na() %>%
      #     .[.$gene_symbol!="",-2]
      #   m_db_3 <-BRCA
      # 
      #   em <- GSEA(score, TERM2GENE=m_db_3, pvalueCutoff=1, pAdjustMethod="BH", eps=0, minGSSize = 5)
      #   write.csv(em, paste0("./",Sys.Date(),"_GSEA_",category,"_",comparison,".csv"))
      #   saveRDS(em, paste0("./",Sys.Date(),"_GSEA_",category,"_",comparison,".RDS"))
      #   em.sig<-em@result[em@result$p.adjust<=1,]
      #   em.sig.plot<-em.sig$ID[em.sig$p.adjust<=1]
      #   for (i in 1:length(em.sig.plot)){
      #     geneSetID<-em.sig.plot[i]
      #     NES<-em@result$NES[em@result$ID==geneSetID]
      #     NES<-format(NES, digits=4)
      #     pvalue<-em@result$p.adjust[em@result$ID==geneSetID]
      #     pvalue<-format(pvalue, digits=4)
      #     FDR<-em@result$qvalues[em@result$ID==geneSetID]
      #     FDR<-format(FDR, digits=4)
      #     p1<-gseaplot2(em, geneSetID=geneSetID, title=geneSetID)+
      #       theme(plot.title = element_text(size=0.5))+
      #       annotate("text", x=0.95, y=0.9, label= paste("NES=",NES), hjust = 1)+
      #       annotate("text", x=0.95, y=0.85, label= paste("adj-p-value=",pvalue), hjust = 1)+
      #       annotate("text", x=0.95, y=0.8, label= paste("FDR=",FDR), hjust = 1)
      #     p1
      #     ggsave(paste0(date,"_GSEA_",category,"_",comparison,"-",str_sub(geneSetID, 10, str_count(geneSetID)),".png"),dpi = 300)
      #   }
      # }
    }

  }
}


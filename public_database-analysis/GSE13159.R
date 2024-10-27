
library(openxlsx)
library("GEOquery")
library(tidyverse)


setwd(r'(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\AML-public_database\GSE13159)')
list.files()

cohort<-"GSE13159"

# read matrix table
filename<-"GSE13159_series_matrix.txt"
gse <- getGEO(filename=filename)

feature<-gse@featureData@data


# filter out the probes of GOI
gene.name<-c("HDAC8")

feature.index<-which(feature$`Gene Symbol` %in% gene.name)
feature.filtered<- feature[feature.index,]
print(feature.filtered$ID)

matrix.filtered<-gse@assayData$exprs[feature.index,] %>% 
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "rowname")


meta<-gse@phenoData@data %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "rowname")

df <-meta %>%
  left_join(., matrix.filtered,  by=c("rowname" = "rowname")) 

unique(df$leukemia_class)

colnames(df)[colnames(df)=="leukemia class:ch1"] <-"leukemia_class"
colnames(df)[colnames(df)=="sample type:ch1"] <-"sample_type"
df$leukemia_class<-factor(df$leukemia_class, levels = c(
  "Non-leukemia and healthy bone marrow",
  "AML with inv(16)/t(16;16)",
  "AML with t(8;21)",
  "AML with t(11q23)/MLL",
  "AML with t(15;17)",
  "ALL with t(12;21)",
  "ALL with t(1;19)",
  "AML with normal karyotype + other abnormalities",
  "ALL with hyperdiploid karyotype",
  "AML complex aberrant karyotype",
  "MDS",
  "T-ALL",
  "c-ALL/Pre-B-ALL with t(9;22)",
  "c-ALL/Pre-B-ALL without t(9;22)",
  "Pro-B-ALL with t(11q23)/MLL",
  "mature B-ALL with t(8;14)",
  "CML",
  "CLL"))

file.name<-paste(Sys.Date(),cohort,gene.name,".xlsx", sep = "-")
print(file.name)
openxlsx::write.xlsx(df,file =file.name)


print(feature.filtered$ID)
# check HDAC8 in different variable
filename<-paste(Sys.Date(),cohort,"box plot-leukemia_class-223345_at(HDAC8)",".png",sep="-")
print(filename)
png(filename, res=300, units="in", height=8, width=10)
df%>%
  ggplot(., aes(leukemia_class, `223345_at`, fill=leukemia_class, color=leukemia_class)) +
  geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.margin = margin(1, 1, 1, 5, "cm"))
graphics.off()  

filename<-paste(Sys.Date(),cohort,"box plot-leukemia_class-223908_at(HDAC8)",".png",sep="-")
print(filename)
png(filename, res=300, units="in", height=8, width=10)
df%>%
  ggplot(., aes(leukemia_class, `223908_at`, fill=leukemia_class, color=leukemia_class)) +
  geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.margin = margin(1, 1, 1, 5, "cm"))
graphics.off() 

filename<-paste(Sys.Date(),cohort,"box plot-leukemia_class-223909_s_at(HDAC8)",".png",sep="-")
print(filename)
png(filename, res=300, units="in", height=8, width=10)
df%>%
  ggplot(., aes(leukemia_class, `223909_s_at`, fill=leukemia_class, color=leukemia_class)) +
  geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=1))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), alpha = 1/5)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.margin = margin(1, 1, 1, 5, "cm"))
graphics.off() 

library(ggplot2)
library(dplyr)

# code 1: basic stats.

# ==================================================================================

# part1: make the summary plot from the home-summarized table 
# read the saved RDS from here:
project_name<-"20230906_HDAC8-134"
# define rMATS output folder:
root_folder<-r'()'
setwd(paste0(root_folder,"\\",project_name))
list.files()

summary_table<-readRDS("2024-01-24-splicing_summary-all.RDS")

ggplot(summary_table, aes(x=FDR, color=cohort)) + 
  geom_density()+
  ggtitle("unpaired")

mat_long<-summary_table %>%
  group_by(cohort, AS_type,  .drop = FALSE) %>%
  filter(PValue <= 0.05) %>%
  filter(FDR <= 0.05) %>%
  count()

mat_long$cohort<-factor(mat_long$cohort, level=c("Ctrl-IR vs Ctrl-NIR","CM-IR vs CM-NIR"))

# dodge bar graph
ggplot(mat_long, aes(x=AS_type, y=n, fill=cohort, label=n)) +
  geom_bar(stat = "identity", position = "dodge")+
  geom_text(size = 3, position = position_dodge(width = 1), vjust=-1)+
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  ggtitle("unpaired: p-value<=0.05, FDR<=0.05")

# stacked bar graph
ggplot(mat_long, aes(x=cohort, y=n, fill=AS_type, label=n)) +
  geom_bar(stat = "identity", position = "stack")+
  geom_text(size = 3, position = position_stack(vjust = 0.5))



#make the summary plot from the summary table from rMATS

project_name<-"20240121_paired"
setwd(paste0("C:\\Users\\yufu\\OneDrive - City of Hope National Medical Center\\rotation\\1_Dr. Kuo\\U2AF1\\RNAseq\\splicing\\rMATS\\",project_name))
list.files()
summary_1<-read.table(".\\1\\final_output\\summary.txt", sep = "\t", header = T, quote = "")
summary_2<-read.table(".\\2\\final_output\\summary.txt", sep = "\t", header = T, quote = "")

mat<-as.data.frame(matrix(nrow = 5, ncol = 0))
mat$type<-summary_1$EventType
mat$WT_IR_vs_NIR<-summary_1$SignificantEventsJCEC
mat$CM_IR_vs_NIR<-summary_2$SignificantEventsJCEC
mat

mat_long <- mat %>% 
  pivot_longer(
    cols = 2:3, 
    names_to = "sample",
    values_to = "events"
  )
unique(mat_long$sample)
mat_long$sample<-factor(mat_long$sample, level=c("WT_IR_vs_NIR","CM_IR_vs_NIR"))


# dodge bar graph
ggplot(mat_long, aes(x=type, y=events, fill=sample, label=events)) +
  geom_bar(stat = "identity", position = "dodge")+
  geom_text(size = 3, position = position_dodge(width = 1), vjust=-1)

# stacked bar graph
ggplot(mat_long, aes(x=sample, y=events, fill=type, label=events)) +
  geom_bar(stat = "identity", position = "stack")+
  geom_text(size = 3, position = position_stack(vjust = 0.5))



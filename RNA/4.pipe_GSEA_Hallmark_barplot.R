#### Bar plot of GSEA result ( NSM vs. Normal, nonNSM vs. Normal, NSM vs. nonNSM) ####

## Import Libary
library(ggplot2)

## Parameters
setwd("/home/jen96/ESCC/RNA/GSEA")
input_path <- getwd()

## Color: 
# NSM  #8ab594
# nonNSM #ada7a0
# Normal  #dbd7c5

# 1. NSM vs Normal
output_filename<-"norm_vst.matrix.filtered.GSEA.NSM_Normal.barplot.pdf"
df<-read.table(paste0(input_path,"/norm_vst.matrix.filtered.GSEA.NSM_normal.filter.txt"),header=T,sep='\t')
color <- ifelse(df$NES > 0, "#8ab594","#dbd7c5")

# 2. nonNSM vs Normal
output_filename<-"norm_vst.matrix.filtered.GSEA.nonNSM_Normal.barplot.pdf"
df<-read.table(paste0(input_path,"/norm_vst.matrix.filtered.GSEA.nonNSM_normal.filter.txt"),header=T,sep='\t')
color <- ifelse(df$NES > 0, "#ada7a0","#dbd7c5")

# 3. NSM vs nonNSM
output_filename<-"norm_vst.matrix.filtered.GSEA.NSM_nonNSM.barplot.pdf"
df<-read.table(paste0(input_path,"/norm_vst.matrix.filtered.GSEA.NSM_nonNSM.filter.txt"),header=T,sep='\t')
color <- ifelse(df$NES > 0, "#8ab594","#ada7a0")

# Saved as Fig. 4B
ggplot(df ,aes(x=reorder(HALLMARK,NES), y=NES))+
  geom_bar(stat='identity',show.legend = F, 
           fill=color,
           color='white')+
  coord_flip()+
  ylab("Normalized Enrichment Score (NES)")+
  xlab("Hallmark")+
  theme_minimal()+
  geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
  geom_text(aes(label = HALLMARK,
                hjust = ifelse(NES < 0, -0.1,1.05),
                vjust = 0.5), size = 4) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank())


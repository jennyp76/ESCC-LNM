#### Unsupervised hierarchical clustering using gene expression profiles of 10 EMT marker genes ####

## Import Library
library(ComplexHeatmap)

## Parameters
setwd("/home/jen96/ESCC/RNA/EMT")
input_path <- getwd()
EMT_marker<-c("VIM","FN1","ZEB1","CDH2","MMP2","SNAI2","ZEB2","SPARC","SNAI1","CDH1")

# Subset expression corresponding to EMT_marker
exp <- read.table(paste0(input_path,"/norm_vst.matrix.filtered.NSM_nonNSM.txt"), header=T, row.names=1 )
df <- data.frame()
for (i in rownames(exp)){
  if( i %in% EMT_marker){
    line <- exp[i,]
    df <- rbind(df,line)}}
names(df)<-names(exp)

# Plotting 
  # Saved as Fig. 4C
data_scaled <- t(scale(t(df)))
column_anno <- HeatmapAnnotation(foo=anno_block(gp=gpar(fill=c("#8ab594","#dbd7c5","#ada7a0")),labels=c("NSM","Normal","nonNSM")))
row_anno <- rowAnnotation(foo=anno_block(gp=gpar(fill=c("#C1ECE4","#FFD0D0"))))
Heatmap(data_scaled,
        clustering_method_rows = "complete",
        column_title = NULL,
        row_title =NULL,
        show_row_names         = T,
        show_column_names      = T,
        cluster_columns        = T,
        cluster_rows           = T,
        top_annotation = column_anno,
        column_km=3,
        right_annotation= row_anno,
        row_km=2)


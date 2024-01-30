#### Quality Control for all transcriptomic data using GEPIA database ####

## Import Library 
library("ComplexHeatmap")
library("DESeq2")
library("dplyr")

## Parameters
setwd("/home/jen96/ESCC/RNA/GEPIA")
input_path <- getwd()
output_path <- getwd()

## Extract VST 
# Read raw count data
countdata <- read.table(paste0(input_path,'/Total_Count.RGadded.marked.FeatureCounts_reverse_min10.matrix.txt'), header = TRUE, fill = TRUE, row.names = 1)

# Generate coldata 
subject <- factor(c(rep("1",2),rep("2",2),rep("3",2),rep("4",2),rep("5",2),rep("6",2),rep("7",2)))
cnd <- factor(c("Normal","Primary","Normal","Primary","Normal","Primary","Normal","Primary","Normal","Primary","Normal","Primary","Normal","Primary"))

coldata <- data.frame(row.names = colnames(countdata),cnd,subject)
coldata$cnd <- relevel(coldata$cnd, ref = "Normal")

# Generate DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ subject + cnd )

# Pre-filtering low count genes
keep <- rowSums(counts(dds)) >= length(colnames(countdata))
dds <- dds[keep,]

# Running DESeq
df.deseq <- DESeq(dds)
resultsNames(df.deseq)

# Extract vst fille 
vst.deseq <- vst(df.deseq, blind=FALSE)
norm_vst <- assay(vst.deseq)
write.table(norm_vst, file = paste0(output_path,'/Total.norm_vst.matrix.txt'),
            append = FALSE, quote = FALSE, sep = "\t", eol = '\n', na = 'NA', row.names = TRUE, col.names = NA )

## Conduct GEPIA 
GEPIA_db <- read.table("./GEPIA_ESCA_lfc_1_0.05_simple.txt", header=T)

data <- read.table(paste0(input_path,"/Total.norm_vst.matrix.GEPIA2.txt"), header=T, row.names=1 ) 
data_scaled <- t(scale(t(data)))
data_scaled <-data_scaled[order(row.names(data_scaled)),]

df = data.frame()
for (i in rownames(data_scaled)){
  if (i %in% GEPIA_df$Gene_Symbol){
    line <-c(i,subset(GEPIA_db,Gene_Symbol == i)$Upregulated)
    df = rbind(df,line)
  }
}

colnames(df) <- c("Gene","Upregulated")
rownames(df) <- df$Gene

df2 <- subset(df, select =-c(Gene))
df2$Upregulated <- factor(df2$Upregulated, levels= c("Normal", "Tumor"))
df2<-as.matrix(df2)

# Saved as Supplementary Fig. S1
ht<-Heatmap(df2, name="Gene_Upregulated",show_row_names= F, col=c("yellow","green"),show_column_names = FALSE,
            heatmap_legend_param = list(title = "Genes Upregulated", title_gp=gpar(fontsize =9, fontface = "bold"), labels_gp=gpar(fontsize=9)))
Heatmap(data_scaled,
        clustering_method_rows = "complete",
        show_row_names         = F,
        cluster_columns        = T,
        cluster_rows           = F,
        row_km = 2,
        column_km = 2,
        column_names_gp = gpar(fontsize = 7),
        heatmap_legend_param = list(title = " Gene Expression", title_gp=gpar(fontsize = 9, fontface = "bold"),labels_gp=gpar(fontsize=9))) + ht
 

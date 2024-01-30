#### Conduct DESeq2 & extract variance stabilizing transformations (VST) file #### 

## Import Library
library("DESeq2")
library("limma")
library("dplyr")
library("vidger")

## Parameters
setwd("/home/jen96/ESCC/RNA/DESeq2")
input_path <- getwd()
output_path <- getwd()

# Read raw count data
countdata <- read.table(paste0(input_path,'/Count.RGadded.marked.FeatureCounts_reverse_min10.matrix.txt'), header = TRUE, fill = TRUE, row.names = 1)

# Generate coldata 
subject <- factor(c(rep("1",1),rep("2",1),rep("3",2),rep("4",2),rep("5",2),rep("6",2),rep("7",1)))
cnd <- factor(c("Normal","non_NSM","Normal","non_NSM","Normal","non_NSM","Normal","NSM","Normal","non_NSM","NSM"))
 
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
write.table(norm_vst, file = paste0(output_path,'/norm_vst.matrix.filtered.NSM_nonNSM.txt'),
            append = FALSE, quote = FALSE, sep = "\t", eol = '\n', na = 'NA', row.names = TRUE, col.names = NA )

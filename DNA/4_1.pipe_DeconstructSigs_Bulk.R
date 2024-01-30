#### Execute deconstructSigs for Total mutations, Pre-metastasis mutations and Post-metastasis mutations detected across all patients####

## Import Libaries
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(graphics)

## Parameters
setwd("/home/jen96/ESCC/DNA/deconstructSigs/")
raw_dir<-getwd()
out_dir <- getwd()


# Step1: Total mutations
Total <- 'Total.RGadded.marked.fixed.merge.rescue.total.deconstructSigs.txt'
sample <- strsplit( Total , split='.', fixed=T)[[1]][1]
assign(Total, read.table(paste0(raw_dir,'/', Total ),header=T))
sigs.input <- mut.to.sigs.input(mut.ref = get(Total) , sample.id = "Sample", chr = "chr", pos = "pos",  ref = "ref", alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
RESULT <- whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.exome.cosmic.v3.may2019,
                          sample.id = sample , contexts.needed = TRUE, tri.counts.method = 'default')
  # Saved as Fig. 3A (top)
deconstructSigs::plotSignatures( RESULT, sig.type = "SBS", sub = "")
  # Saved as Supplementary Fig. S3 (right)
deconstructSigs::makePie(RESULT,sub = '',v3=T)


# Step2: Pre-metastasis mutations  
Pre_met <- 'Pre.RGadded.marked.fixed.merge.rescue.pre.deconstructSigs.txt'
sample <- strsplit(Pre_met , split='.', fixed=T)[[1]][1]
assign(Pre_met, read.table(paste0(raw_dir,'/', Pre_met ),header=T))
sigs.input <- mut.to.sigs.input(mut.ref = get(Pre_met) , sample.id = "Sample", chr = "chr", pos = "pos",  ref = "ref", alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
RESULT <- whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.exome.cosmic.v3.may2019,
                          sample.id = sample , contexts.needed = TRUE, tri.counts.method = 'default')

  # Saved as Fig. 3A (bottom, left)
deconstructSigs::plotSignatures( RESULT, sig.type = "SBS", sub = "")
  # Saved as Supplementary Fig. S3 (middle) 
deconstructSigs::makePie(RESULT,sub = '',v3=T)


# Step3: Post-metastasis mutations
Post_met <- 'Post.RGadded.marked.fixed.merge.rescue.post.deconstructSigs.txt'
sample <- strsplit( Post_met , split='.', fixed=T)[[1]][1]
assign(Post_met, read.table(paste0(raw_dir,'/', Post_met ),header=T))
sigs.input <- mut.to.sigs.input(mut.ref = get(Post_met) , sample.id = "Sample", chr = "chr", pos = "pos",  ref = "ref", alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
RESULT <- whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.exome.cosmic.v3.may2019,
                          sample.id = sample , contexts.needed = TRUE, tri.counts.method = 'default')

  # Saved as Fig. 3A (bottom, right)
deconstructSigs::plotSignatures( RESULT, sig.type = "SBS", sub = "")
  # Saved as Supplementary Fig. S3 (left)
deconstructSigs::makePie(RESULT,sub = '',v3=T)


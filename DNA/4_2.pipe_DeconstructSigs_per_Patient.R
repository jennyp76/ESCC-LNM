#### Execute deconstructSigs for Pre-metastasis mutations and Post-metastasis mutations per Patient ####

## Import Libaries
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(graphics)
library(ggplot2)
library(ggpubr)
library(dplyr)

## Parameters
setwd("/home/jen96/ESCC/DNA/deconstructSigs/")
raw_dir<-getwd()
out_dir <- getwd()
wilcox_input <- data.frame()
wilcox_input <- rbind(wilcox_input, list("SBS_signature","Location","Proportion","Patient"))
colnames(wilcox_input)<-c("SBS_signature","Location","Proportion","Patient")

# Step1. Pre-metastasis mutations per Patient
txt_list <- list.files(raw_dir, pattern = '_pre.RGadded.marked.fixed.merge.rescue.pre.deconstructSigs.txt')
for (idx in 1:length(txt_list)) {
  sample <- strsplit(txt_list[idx], split='.', fixed=T)[[1]][1]
  patient <-strsplit(sample, split='_', fixed = T)[[1]][1]
  assign(txt_list[idx], read.table(paste0(raw_dir,'/', txt_list[idx]),header=T))
  sigs.input <- mut.to.sigs.input(mut.ref = get(txt_list[idx]), sample.id = "Sample", chr = "chr", pos = "pos", ref = "ref", alt = "alt",
                                  bsg = BSgenome.Hsapiens.UCSC.hg38)
  RESULT <- whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.exome.cosmic.v3.may2019,
                            sample.id = sample , contexts.needed = TRUE, tri.counts.method = 'default')

  wilcox_input <- rbind(wilcox_input, list("SBS1","Pre",RESULT[['weights']]['SBS1'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS2","Pre",RESULT[['weights']]['SBS2'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS5","Pre",RESULT[['weights']]['SBS5'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS6","Pre",RESULT[['weights']]['SBS6'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS13","Pre",RESULT[['weights']]['SBS13'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS16","Pre",RESULT[['weights']]['SBS16'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS18","Pre",RESULT[['weights']]['SBS18'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS25","Pre",RESULT[['weights']]['SBS25'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS29","Pre",RESULT[['weights']]['SBS29'][[1]][1],patient))
}

# Step2. Post-metastasis mutations per Patient 
txt_list <- list.files(raw_dir, pattern = '_post.RGadded.marked.fixed.merge.rescue.post.deconstructSigs.txt')
for (idx in 1:length(txt_list)) {
  sample <- strsplit(txt_list[idx], split='.', fixed=T)[[1]][1]
  patient <- strsplit(sample, split='_', fixed = T)[[1]][1]
  assign(txt_list[idx], read.table(paste0(raw_dir,'/', txt_list[idx] ),header=T))
  sigs.input <- mut.to.sigs.input(mut.ref = get(txt_list[idx]), sample.id = "Sample", chr = "chr", pos = "pos",  ref = "ref", alt = "alt",
                                  bsg = BSgenome.Hsapiens.UCSC.hg38)
  RESULT <- whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.exome.cosmic.v3.may2019,
                            sample.id = sample , contexts.needed = TRUE, tri.counts.method = 'default')
    
  wilcox_input <- rbind(wilcox_input, list("SBS1","Post",RESULT[['weights']]['SBS1'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS2","Post",RESULT[['weights']]['SBS2'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS5","Post",RESULT[['weights']]['SBS5'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS6","Post",RESULT[['weights']]['SBS6'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS13","Post",RESULT[['weights']]['SBS13'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS16","Post",RESULT[['weights']]['SBS16'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS18","Post",RESULT[['weights']]['SBS18'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS25","Post",RESULT[['weights']]['SBS25'][[1]][1],patient))
  wilcox_input <- rbind(wilcox_input, list("SBS29","Post",RESULT[['weights']]['SBS29'][[1]][1],patient))
}

# Step3. Wilcox test - Pre-metastasis mutations vs. Post-metastasis mutations
wilcox_input <- wilcox_input[-1,]
wilcox_input$Proportion<-as.numeric(wilcox_input$Proportion)

data <-wilcox_input
data$Location <- factor(data$Location, levels= c('Pre','Post'))

  # Saved as Proportion_of_SBS_signature.pdf
ggplot(data, aes(x=SBS_signature, y = Proportion,fill=Location))+
  geom_boxplot(position=position_dodge(0.9))+
  geom_point(aes(fill=Location,group=Location),shape=10, position = position_jitterdodge(jitter.width = .4, dodge.width = 1))+
  scale_shape_manual(values =c(20,23))+
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size=10), axis.title.y = element_text(size=12),axis.text.y = element_text(size=10))+
  stat_compare_means(aes(group = Location))+
  scale_fill_brewer(palette = "Accent")

# Step4. Subset only significant signatures
subset_data <- data %>%
  filter(SBS_signature=='SBS16' |SBS_signature=='SBS6' |SBS_signature=='SBS5'| SBS_signature=='SBS18' )
subset_data$SBS_signature <- factor(subset_data$SBS_signature, levels= c('SBS16','SBS6','SBS5','SBS18'))
subset_data$Location <- factor(subset_data$Location, levels= c('Pre','Post'))

  # Saved as Fig. 3B
ggplot(subset_data, aes(x=SBS_signature, y = Proportion,fill=Location))+
  scale_y_continuous(name = "Proportion of SBS signature")+
  geom_boxplot(outlier.colour="red", outlier.shape=NA,outlier.size=4,position=position_dodge(0.6),alpha=0.85,width=0.5,lwd=0.2,fatten=3)+
  geom_point(aes(fill=Location,group=Location ),size=0.9,shape=21, position = position_jitterdodge(jitter.width = 0.3, dodge.width =0.6 ))+
  stat_compare_means(aes(group = Location), label ="p.format")+
  scale_fill_manual(values = c("Pre"="#143d59", "Post"="#f4b41a"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


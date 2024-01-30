#### Calculate Spearman Correlation Coefficient between ToM and MARM ####

## Import Libaries
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(broom)

## Parameters
setwd("/home/jen96/ESCC/DNA/ToM_MARM")
setwd("/home/jen96/ESCC/DNA")
output_path <- getwd()

data <-as.data.frame(read_table("Calculation_of_ToM_MARM.txt"))
input <- data[,c(1,4,7)]
colnames(input)<-c('Sample','ToM','MARM')
input[,'Patient']<-rep(c('P1','P2','P3','P4','P5','P6','P7','P8','P9','P10'),c(2,1,2,3,3,3,3,2,3,3))                    
input$Patient <- factor(input$Patient, levels =c("P1","P2","P3",'P4','P5',"P6","P7",'P8','P9','P10'))

# Step1: Calculate Spearman correlation coefficient 
cor.test(input$ToM , input$MARM, method = "spearman", exact=FALSE)

# Step2: Scatter Plot with regression
  # Saved as Fig. 3C 
ggplot(input, mapping=aes(x=ToM,y=MARM))+
  geom_point(aes(color=factor(Patient)),size = 3,shape=18) +
  geom_smooth(method=lm, linewidth=0.5)+
  theme_bw() +
  theme(axis.line = element_line(color='black'),plot.background = element_blank(),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  labs(x= 'Timing of metastastsis(ToM)',y='Mutation accumulation rate after metastasis(MARM)', color="Patient")+
  scale_color_brewer( palette = "Paired")+
  scale_x_continuous(breaks=c(0.3,0.5,0.7))+
  scale_y_continuous(breaks=seq(-3,15,3))

# Step3: Scatter Plot with line among patients
  # Saved as Supplementary Fig. S4
ggplot(input, mapping=aes(x=ToM,y=MARM))+
  geom_point(aes(color=factor(Patient)),size = 3,shape=18) +
  geom_line(aes(group=Patient,color=factor(Patient)),linewidth=0.7)+
  theme_bw() +
  theme(axis.line = element_line(color='black'),plot.background = element_blank(),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  labs(x= 'Timing of metastastsis(ToM)',y='Mutation accumulation rate after metastasis(MARM)', color="Patient")+
  scale_color_brewer( palette = "Paired")+
  scale_x_continuous(breaks=c(0.25,0.5,0.75))+
  scale_y_continuous(breaks=seq(0,15,1.5))





  

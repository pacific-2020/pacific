#!/usr/bin/env R

#This script is designed to create a bubble plot of the proportion of predicted reads >95% cutoff 
#per class from PACIFIC summary.txt files. It also plots the error profile vs total number of reads
# of negative controls using this data.

#Load required libraries
library(tidyverse)

#Load dataframe
dat <- read_tsv("pacific/figures/rnaseqtests.txt", col_names=T)
colnames(dat) <- c("filename", "class", "predicted_reads", "perc_predicted_reads", "predicted_reads_95","perc_predicted_reads_95")

#Create cutoff data frame. Cutoffs were derived from FPR tests
co <- tibble(class = unique(dat$class), cutoff = c(0.0004,0.0506,0.0014,0.0002,0.0466,NA,NA,NA)) %>% 
  filter(!class %in% c("Human", "Discarded", "rc_discarded"))

#Join cutoff information
df <- full_join(dat,co, by="class")

#Create column stating whether files are positive for a virus
df$classification <- ifelse(df$perc_predicted_reads_95>df$cutoff, "Positive", "Negative")

#Remove text
df$filename <- gsub("_.*","",df$filename)
df$filename <- gsub("\\.fa.*","",df$filename)

##############
#Plots for positive and negative controls

negmd <- read_csv("pacific/metadata/rnaseqtests.negsamples.SraRunTable.20200624.csv", col_names=T)
posmd <- read_csv("pacific/metadata/rnaseqtests.possamples.SraRunTable.20200624.csv", col_names=T)

#Filter and label for negative and positive samples
neg <- semi_join(df,negmd,by=c("filename"="Run")) %>%
  mutate(filename=paste(filename,"neg",sep="."))
  #mutate(type="SARS-CoV-2 -")

pos <- semi_join(df,posmd,by=c("filename"="Run")) %>%
  mutate(filename=paste(filename,"pos",sep="."))
  #mutate(type="SARS-CoV-2 +")
  
#Join negative and positive tables
negpos <- rbind(neg,pos)

#Create bubble plot of classifications for positive and negative controls
negpos %>% 
  select(filename,class,perc_predicted_reads_95, classification,cutoff) %>%
  filter(!class %in% c("Human", "Discarded", "rc_discarded")) %>%
  arrange(desc(classification,filename)) %>%
  mutate(filename = fct_rev(as_factor(filename))) %>%
  ggplot() +
  geom_point(aes(x=class,y=filename, size=cutoff), color="black", alpha=0.3) +
  geom_point(aes(x=class,y=filename, size=perc_predicted_reads_95, color=classification), alpha=1) +
  scale_color_manual(name="Classification", values=c("dodgerblue", "black"), breaks=c("Positive","Negative"))+
  labs(y = "Sample", x="Virus class")+
  guides(size=guide_legend(title="PR (%)"))+
  theme_bw()+
  theme(text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1))
  #scale_x_discrete(breaks=c("DRR028492", "SRR10303453", "SRR4895842", "SRR5375307", 
  #  "SRR5515378", "SRR6116514", "SRR6411409", "SRR8927151", "SRR8928257", "SRR9274276", 
  #  "SRR10971381", "SRR11412227", "SRR11412228", "SRR11412229","SRR11412230"),
  #  labels=c("SARS-CoV-2 -", "SARS-CoV-2 -", "SARS-CoV-2 -","SARS-CoV-2 -",
  #           "SARS-CoV-2 -","SARS-CoV-2 -","SARS-CoV-2 -","SARS-CoV-2 -",
  #           "SARS-CoV-2 -","SARS-CoV-2 -","SARS-CoV-2 +","SARS-CoV-2 +",
  #           "SARS-CoV-2 +","SARS-CoV-2 +","SARS-CoV-2 +"))
  #theme(axis.title=element_text(size=10), axis.text.x=element_text(angle=90,size=10),axis.text.y=element_text(size=15), legend.text = element_text(size = rel(1.5)), legend.title = element_text(size = rel(1.5)))
  

#Plot error vs total number of reads for negative control samples
#Is there a correlation between negative controls and total reads?
neg %>% 
  group_by(filename) %>% 
  summarise(totalreads= sum(predicted_reads), error=sum(perc_predicted_reads_95)-sum(perc_predicted_reads_95[6:8])) %>% 
  ggplot(aes(totalreads,error)) +
  geom_point(alpha=0.5)+
  labs(y = "Error", x="Total number of reads in sample")

###########
#Plots for real dataset experiments (Wesolowska-Andersen et al. 2017)

#Load manuscript metadata table
md <- read_tsv("pacific/metadata/Wesolowska-Andersen2017.combined.tsv", col_names=T)

#Join md and df, retaining sample only in md
combdf <- left_join(md,df,by=c("run"="filename"))
wcombdf <- combdf %>% 
  select(run,class, sampleID, classification,perc_predicted_reads_95,cutoff,qpcr_species,rnaseq_species,rnaseq_scaledcoverage) %>%
  filter(!class %in% c("Human", "Discarded", "rc_discarded")) %>%
  filter(classification %in% "Positive")

#Extract positive samples
#md %>% filter(run %in% c(wcombdf %>% pull(run))) %>%
#  select(run,sampleID)

#Create bubble plot of classifications for md samples
combdf %>% 
  select(run,class,perc_predicted_reads_95, classification,cutoff) %>%
  filter(!class %in% c("Human", "Discarded", "rc_discarded")) %>%
  drop_na %>%
  arrange(desc(classification,run)) %>%
  mutate(run = fct_rev(as_factor(run))) %>%
  ggplot() +
  geom_point(aes(class,run, size=cutoff), color="black", alpha=0.3) +
  geom_point(aes(class,run, size=perc_predicted_reads_95, color=classification), alpha=1) +
  scale_color_manual(name="Classification", values=c("dodgerblue", "black"), breaks=c("Positive","Negative"))+
  labs(y = "Sample", x="Class")+
  guides(size=guide_legend(title="PR (%)"))+
  theme_bw()+
  theme(text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1))

##############
#Notes

#Create heatmap of samples v classifcation
#df %>% 
#  select(filename,class,perc_predicted_reads_95, classification,cutoff) %>%
#  filter(!class %in% c("Human", "Discarded", "rc_discarded")) %>%
#  ggplot() +
#  geom_tile(aes(filename,class,fill=perc_predicted_reads_95/cutoff)) +
#  theme(axis.text.x=element_text(angle=90))+
#  scale_fill_continuous(name = "PR>0.95/Cutoff")+
#  labs(y = "Class", x="Sample")

#Create boxplot of RNA seq tests
#df %>% 
#  select(filename,class,perc_predicted_reads_95, type) %>%
#  filter(!class %in% c("Human", "Discarded", "rc_discarded")) %>%
#  ggplot(aes(fill=type)) +
#  geom_boxplot(aes(class,perc_predicted_reads_95), alpha=1) + 
#  #geom_point(aes(class,perc_predicted_reads_95), alpha=0.5) +
#  geom_hline(data = co, aes(yintercept=cutoff), color="red") +
#  theme(legend.position = "bottom") +
#  labs(y = "Predicted reads >95% cutoff (%)", x="Class")+ 
#  facet_wrap(.~class, scales="free")

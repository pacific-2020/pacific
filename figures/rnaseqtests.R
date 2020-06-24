#!/usr/bin/env R

#This script is designed to create a bubble plot of the proportion of predicted reads >95% cutoff 
#per class from PACIFIC summary.txt files. It also plots the error profile vs total number of reads
# of negative controls using this data.

#Load required libraries
library(tidyverse)

#Load dataframe
dat <- read_tsv("pacific/figures/rnaseq.tests.txt", col_names=T)
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
df$filename <- gsub("\\.fasta.gz","",df$filename)

#Create bubble plot of classifications
df %>% 
  select(filename,class,perc_predicted_reads_95, classification,cutoff) %>%
  filter(!class %in% c("Human", "Discarded", "rc_discarded")) %>%
  ggplot() +
  geom_point(aes(class,filename, size=cutoff), color="black", alpha=0.3) +
  geom_point(aes(class,filename, size=perc_predicted_reads_95, color=classification), alpha=1) +
  scale_color_manual(name="Classification", values=c("dodgerblue", "black"), breaks=c("Positive","Negative"))+
  labs(y = "Sample", x="Class")+
  guides(size=guide_legend(title="PR (%)"))+
  theme_bw()

#Plot error vs total number of reads for negative control samples
df %>% 
  group_by(filename) %>% 
  filter(!classification %in% "Positive") %>%
  summarise(totalreads= sum(predicted_reads), error=sum(perc_predicted_reads_95)-sum(perc_predicted_reads_95[6:8])) %>% 
  ggplot(aes(totalreads,error)) +
  geom_point(alpha=0.5)+
  labs(y = "Error", x="Total number of reads in sample")

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

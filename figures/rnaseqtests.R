#!/usr/bin/env R

#This script is designed to create boxplots of the proportion of predicted reads >95% cutoff 
#per class from PACIFIC summary.txt files. It also plots the error profile vs total number of reads
#using this data.

#Load required libraries
library(tidyverse)

#Load dataframe
df <- read_tsv("pacific/figures/rnaseq.tests.txt", col_names=T)
colnames(df) <- c("filename", "class", "predicted_reads", "perc_predicted_reads", "predicted_reads_95","perc_predicted_reads_95")

#Create cutoff data frame. Cutoffs were derived from FPR tests
co <- tibble(class = unique(df$class), cutoff = c(0.0004,0.0506,0.0014,0.0002,0.0466,NA,NA,NA)) %>% 
  filter(!class %in% c("Human", "Discarded", "rc_discarded"))

#Create boxplot of RNA seq tests
df %>% 
  select(filename,class,perc_predicted_reads_95) %>%
  filter(!class %in% c("Human", "Discarded", "rc_discarded")) %>%
  ggplot() +
  geom_boxplot(aes(class,perc_predicted_reads_95), alpha=1) + 
  geom_point(aes(class,perc_predicted_reads_95), alpha=0.5) +
  geom_hline(data = co, aes(yintercept=cutoff, color="red")) +
  theme(legend.position = "none") +
  labs(y = "Predicted reads >95% cutoff (%)", x="Class")+ 
  facet_wrap(.~class, scales="free")

#Plot error vs total number of reads
df %>% 
  group_by(filename) %>% 
  summarise(totalreads= sum(predicted_reads), error=sum(perc_predicted_reads_95)-sum(perc_predicted_reads_95[6:8])) %>% 
  ggplot(aes(totalreads,error)) +
  geom_point(alpha=0.5)+
  labs(y = "Error", x="Total number of reads in sample")

#!/usr/bin/env R

#This script is designed to create bar plots of FNR reads

#Load required libraries
library(tidyverse)

#Load dataframe containing concatenated PACIFIC summaries
dat <- read_delim("pacific/figures/fnrtestresults.txt", delim=":", col_names=F)
colnames(dat) <- c("reads", "set", "class", "type", "FNR")

p1 <- dat %>%
  group_by(class,type,FNR) %>%
  summarise(avg=mean(reads)) %>%
  filter(!FNR %in% c("Discarded")) %>%
  filter(!type %in% "INDEL") %>%
  ggplot(aes(x=type, y=avg, fill=FNR)) +
  geom_bar(stat="identity", position = "stack")+
  labs(x = "True class", y="Mean false negative reads")+
  theme_bw() +
  labs(fill = "Predicted class") +
  facet_grid(~ class) +
  theme(axis.text.x=element_text(angle=30,hjust=1))

p2 <- dat %>%
  group_by(class,type,FNR) %>%
  summarise(avg=mean(reads)) %>%
  filter(!FNR %in% c("Discarded", "rc_discarded", "Human")) %>%
  ggplot(aes(x=class, y=avg, fill=FNR)) +
  geom_bar(stat="identity", position = "stack") +
  labs(x = "True class", y="Mean false negative reads")+
  theme_bw() +
  labs(fill = "Predicted class")+
  theme(axis.text.x=element_text(angle=30,hjust=1))

grid.arrange(p1,p2, ncol=2)

####
#Notes

#Filtered for >=0.95 cutoff
#dat <- read_delim("pacific/figures/fnrtestresults.0.95.txt", delim=":", col_names=F)

#Total reads
#p1 <- dat %>%
#  filter(!FNR %in% c("Discarded")) %>%
#  ggplot(aes(x=class, y=reads, fill=FNR)) +
#  geom_bar(stat="identity", position = "stack")+
#  labs(x = "True class", y="No. of reads")+
#  theme_bw() +
#  labs(fill = "Predicted class")

#p2 <- dat %>%
#  filter(!FNR %in% c("Discarded", "rc_discarded", "Human")) %>%
#  ggplot(aes(x=class, y=reads, fill=FNR)) +
#  geom_bar(stat="identity", position = "stack") +
#  labs(x = "True class", y="No. of reads")+
#  theme_bw() +
#  labs(fill = "Predicted class")+
#  theme(axis.text.x=element_text(angle=30,hjust=1))

#p3 <- dat %>%
  #filter(!FNR %in% c("Discarded")) %>%
  #ggplot(aes(x=class, y=reads, fill=type)) +
  #geom_bar(stat="identity", position = "stack")+
  #labs(x = "True class", y="No. of reads")+
  #theme_bw() +
  #labs(fill = "Type")

#p4 <- dat %>%
#  filter(!FNR %in% c("Discarded", "rc_discarded", "Human")) %>%
#  ggplot(aes(x=class, y=reads, fill=type)) +
#  geom_bar(stat="identity", position = "stack") +
#  labs(x = "True class", y="No. of reads")+
#  theme_bw() +
#  labs(fill = "Type")

#grid.arrange(p1,p2,p3,p4, nrow=2, ncol=2)
#grid.arrange(p5,p2, ncol=2)
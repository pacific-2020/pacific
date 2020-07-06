#!/usr/bin/env R

#This script is designed to create bar plots of FNR reads

#Load required libraries
library(tidyverse)

#Load dataframe containing concatenated PACIFIC summaries
dat <- read_delim("pacific/figures/fnrtestresults95.txt", delim=":", col_names=F)
colnames(dat) <- c("set", "class", "type", "FNR")

p1 <- dat %>%
  group_by(class,type,FNR) %>%
  summarise(avg=n()/100) %>%
  filter(!FNR %in% c("Discarded")) %>%
  filter(!class %in% c("Human")) %>%
  ggplot(aes(x=type, y=avg, fill=FNR)) +
  geom_bar(stat="identity", position = "stack")+
  labs(x = "True class", y="Mean false negative reads")+
  theme_bw() +
  labs(fill = "Predicted class") +
  facet_grid(~ class) +
  theme(axis.text.x=element_text(angle=30,hjust=1))

p2 <- dat %>%
  group_by(class,type,FNR) %>%
  summarise(avg=n()/100) %>%
  filter(!FNR %in% c("Discarded", "rc_discarded", "Human")) %>%
  filter(!class %in% c("Human")) %>%
  ggplot(aes(x=class, y=avg, fill=FNR)) +
  geom_bar(stat="identity", position = "stack") +
  labs(x = "True class", y="Mean false negative reads")+
  theme_bw() +
  labs(fill = "Predicted class")+
  theme(axis.text.x=element_text(angle=30,hjust=1))

grid.arrange(p1,p2, ncol=2)
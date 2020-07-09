#!/usr/bin/env R

#This script is designed to convert concatenated PACIFIC summaries to metadata

#Load required libraries
library(tidyverse)

#Load dataframe containing concatenated PACIFIC summaries
dat <- read_tsv("pacific/figures/rnaseqtests.txt", col_names=T)
colnames(dat) <- c("filename", "class", "preads", "pprop", "p95reads", "p95prop")

#Separate above and below thresholds
tmp1 <- dat %>%
  select(filename,class, preads, pprop) %>%
  mutate(threshold_0.95="below") %>%
  select(filename,class, threshold_0.95, preads, pprop) %>%
  setNames(c("SRA_ID","class","pp_0.95","reads","prop_totalreads"))

tmp2 <- dat %>%
  select(filename,class, p95reads, p95prop) %>%
  mutate(threshold_0.95="above") %>%
  select(filename,class, threshold_0.95, p95reads, p95prop) %>%
  setNames(c("SRA_ID","class","pp_0.95","reads","prop_totalreads"))

#Concatenate files and remove file extensions
combined.df <- rbind(tmp2,tmp1)
combined.df$SRA_ID <- gsub("\\.fa.*","",combined.df$SRA_ID)
combined.df <- combined.df %>%
  arrange(SRA_ID, class)

write_tsv(combined.df, "pacific/metadata/pacificsummaries.realdatasets.txt")



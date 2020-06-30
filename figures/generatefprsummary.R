#!/usr/bin/env Rscript

library(Biostrings)
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
testclass = args[1]
testdata = args[2]
x <- readDNAStringSet(paste("~/pacificfpr/",testclass,"/pacificoutput_rc.pacificoutput_fprtest", testdata, ".", testclass, ".genomes.fasta.gz", sep=""))
x <- names(x)
x <- separate(data.frame(x = x), x, paste("c",1:11,sep = ""), sep = ":")
x$c8 <- as.numeric(x$c8)
x$c10 <- as.numeric(x$c10)
x <- x %>% mutate(c9 = case_when(c8 != -1 & c8 <= 0.95 ~ "pDiscarded", TRUE ~ c9))
x <- x %>% mutate(c11 = case_when(c10 != -1 & c10 <= 0.95 ~ "pDiscarded", TRUE ~ c11))
x <- x %>% mutate(c11 = case_when(c9 != c11 & c10 > 0.95 & c8 > 0.95 ~ "cDiscarded", TRUE ~ c11))
x <- data.frame(table(x$c1,x$c2,x$c7,x$c9,x$c11)) %>% filter(Freq>0) %>% mutate(testclass=testclass)
write.table(x, file = paste("~/pacificfpr/",testclass,"/summary.fprtest", testdata, ".", testclass, ".txt", sep=""), sep = "\t", col.names=F, row.names=F, quote=F)

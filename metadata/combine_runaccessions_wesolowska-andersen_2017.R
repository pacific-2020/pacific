#Load libraries
library(tidyverse)
library(rvest)

#Import manuscript table into R
url <- "https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1140-8/tables/1"
tab <- read_html(url) %>%
  html_nodes("table") %>%
  html_table(fill=T) %>%
  .[[1]] %>%
  setNames(c("subject","qpcr_species","qpcr_count","qpcr_ccl8","qpcr_cxcl11","rnaseq_totalreads","rnaseq_virusreads","rnaseq_percviralreads","rnaseq_species","rnaseq_strain_serotype","rnaseq_ncbiref","rnaseq_percmap2ncbiref","rnaseq_scaledcoverage"))

#Wrangle table
tab <- tibble(tab[-(1:2), , drop = FALSE])
tab <- tab %>%
  mutate_at(c(5:7,12), parse_number)

#Load NCBI accessions
dat <- read_csv("pacific/metadata/Wesolowska-Andersen2017.SraRunTable.txt", col_names=T)
df <- dat %>% 
  select(Run,BioProject,BioSample,Experiment,`Sample Name`) %>%
  set_names("run","bioProject","bioSample","experiment","sampleID")

#Join accessions and table
out <- full_join(df,tab, by=c("sampleID"="subject"))
write_tsv(out,"pacific/metadata/Wesolowska-Andersen2017.combined.tsv")
#!/usr/bin/env R

#This script is designed to create a bubble plot of the proportion of predicted reads >95% cutoff 
#per class from PACIFIC summary.txt files. It also plots data from Wesolowska-Andersen et al. 2017
#using PACIFIC and BWA-MEM.

#Load required libraries
library(tidyverse)
library(reshape2)
library(gridExtra)
library(cowplot)

#Load dataframe containing concatenated PACIFIC summaries
dat <- read_tsv("pacific/figures/rnaseqtests.txt", col_names=T)
colnames(dat) <- c("filename", "class", "predicted_reads", "perc_predicted_reads", "predicted_reads_95","perc_predicted_reads_95")

#Create cutoff data frame. Cutoffs were derived from FPR tests
#Cutoffs were determined as three standard deviations away from the mean in FPR tests
co <- tibble(class = unique(dat$class), 
             cutoff = c(0.0002540988,0.0427490964,0.0008649834,0.0001261678,0.0452518826,NA,NA,NA)) %>% 
  filter(!class %in% c("Human", "Discarded", "rc_discarded"))

#Join cutoff information
df <- full_join(dat,co, by="class")

#Create column stating whether files are positive for a virus
df$classification <- ifelse(df$perc_predicted_reads_95>df$cutoff, "Positive", "Negative")

#Rename text
df$filename <- gsub("_.*","",df$filename)
df$filename <- gsub("\\.fa.*","",df$filename)
df$class <- gsub("SARS-CoV-2", "Sars_cov_2", df$class)

##############
#Figure 5
#Plots for positive and negative controls

#Load required metadata
negmd <- read_csv("pacific/metadata/rnaseqtests.negsamples.SraRunTable.20200624.csv", col_names=T)
posmd <- read_csv("pacific/metadata/rnaseqtests.possamples.SraRunTable.20200624.csv", col_names=T)

#Filter and label for negative and positive samples
neg <- semi_join(df,negmd,by=c("filename"="acc")) %>%
  mutate(filename=paste(filename,"neg",sep="."))

pos <- semi_join(df,posmd,by=c("filename"="Run")) %>%
  mutate(filename=paste(filename,"pos",sep="."))
  
#Join negative and positive tables
negpos <- rbind(neg,pos)
negpos$class <- gsub("SARS-CoV-2", "Sars_cov_2", negpos$class)

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
  theme(text=element_text(size=16),
        axis.text.x=element_text(angle=30,hjust=1))

###########
#Figure 6 and Supplementary Figure 2
#Plots for real dataset experiments (Wesolowska-Andersen et al. 2017)

#Load manuscript metadata table
md <- read_tsv("pacific/metadata/Wesolowska-Andersen2017.combined.tsv", col_names=T)

#Join md and df, retaining sample only in md
combdf <- left_join(md,df,by=c("run"="filename"))

#Load dataframe of PACIFIC and BWA. Generate overlaps
sam <- read_tsv("pacific/figures/samreadids.unique.txt", col_names=F)
fa <- read_tsv("pacific/figures/pacificfareadids.txt", col_names=F)

colnames(fa) <-  c("run", "readid", "class")
colnames(sam) <-  c("run", "readid", "class")

#Text processing
fa$run <- gsub("pacificoutput_", "", fa$run)
sam$class <- gsub("Cornidovirineae", "Coronaviridae", sam$class)


#Separate by class for correlation test
bwa <- sam %>% 
  group_by(run,class) %>% 
  summarise(bwa=n())

pacific <- fa %>% 
  group_by(run,class) %>% 
  summarise(fa=n())

bwa$class <- gsub("rhinovirus", "Rhinovirus", bwa$class )

pacbwa <- left_join(df,bwa,by=c("filename"="run", "class"="class")) %>% 
  filter(filename %in% md$run) %>%
  select(filename,class,predicted_reads_95, bwa) %>%
  replace(is.na(.), 0) %>%
  filter(!class %in% c("Human", "Discarded", "rc_discarded"))

#Correlation test - Spearman's correlation is used as there are outliers on the data 
#Therefore, Pearson would be inappropriate

result <- NULL
for (c in unique(pacbwa$class))
{
  tmp <- pacbwa %>% 
    filter(class %in% c)
  o <- cor(tmp$predicted_reads_95,tmp$bwa, method="spearman")
  result[c] <- o
}

#As all SARS-CoV-2 is zero, assign correlation of 1
result[1] <- 1
mean(result)

#Report positive classes
bwaj <- left_join(bwa,md, by="run") %>%
  select(run,class,bwa,qpcr_species,rnaseq_species) %>%
  arrange(desc(qpcr_species,rnaseq_species, bwa))
pacificj <- left_join(pacific,md, by="run") %>%
  select(run,class,fa,qpcr_species,rnaseq_species) %>%
  arrange(desc(qpcr_species,rnaseq_species, fa))

#Prepare dataframes
s <- sam %>% 
  mutate(id = paste(run,readid,class, sep=":"), type="sam") %>%
  select(run, id,type)

f <- fa %>% 
  mutate(id = paste(run,readid,class, sep=":"), type="fa") %>%
  select(run,id,type)

j <- full_join(f,s, by=c("run","id"))
colnames(j) <- c("run", "id", "fa.pres", "sam.pres")

#Join dataframes
both <- j %>% 
  drop_na() %>% 
  group_by(run) %>% 
  summarise(overlap=n())

ftotal <- j %>% 
  filter(fa.pres %in% "fa") %>% 
  group_by(run) %>% 
  summarise(pacific=n())

stotal <- j %>% 
  filter(sam.pres %in% "sam") %>% 
  group_by(run) %>% 
  summarise(bwa=n())

tmp <- full_join(stotal,both)

#Create summary table
summary <- full_join(ftotal,tmp) %>%
  replace(is.na(.), 0) %>% 
  mutate(BWA_only = bwa-overlap, PACIFIC_only=pacific-overlap)
  
#Combine dfs
combdf<- full_join(combdf,summary, by="run")

#Convert NAs to zero
wcombdf <- combdf %>%
  select(run,class,perc_predicted_reads_95, classification,cutoff,BWA_only,PACIFIC_only,overlap) %>%
  filter(!class %in% c("Human", "Discarded", "rc_discarded")) 

#Get rid of rows with all NAs
wcombdf <- wcombdf[rowSums(is.na(wcombdf)) != ncol(wcombdf), ]

#Prepare dataframe for plotting
wcombdf <- wcombdf %>%
  mutate(total = BWA_only+PACIFIC_only+overlap)%>%
  drop_na() %>%
  arrange(desc(total,classification)) %>%
  mutate(run = fct_rev(as_factor(run)))

#Filter for only positive samples - Figure 6
pcombdf <- wcombdf %>%
  filter(run %in% wcombdf[wcombdf$classification=="Positive",]$run)

#Create bubble plot of classifications for md samples
p1 <- wcombdf %>% 
  select(run,class,perc_predicted_reads_95, classification,cutoff) %>%
  ggplot() +
  geom_point(aes(class,run, size=cutoff), color="black", alpha=0.3) +
  geom_point(aes(class,run, size=perc_predicted_reads_95,color=classification), alpha=1) +
  scale_color_manual(name="Classification", values=c("dodgerblue", "black"), breaks=c("Positive","Negative"))+
  labs(y = "Sample", x="Virus class")+
  guides(size=guide_legend(title="PR (%)"))+
  theme_bw()+
  theme(
        axis.text.x=element_text(angle=0,hjust=1),
        axis.title=element_blank(),
        legend.position = "bottom",
        )

p2 <-  wcombdf %>%
  select(run,classification,BWA_only,PACIFIC_only,overlap) %>%
  gather(overlap,BWA_only,PACIFIC_only, key="Key", value="value") %>%
  mutate(Key=factor(Key, levels=c("BWA_only", "PACIFIC_only", "overlap"))) %>%
  drop_na() %>% 
  ggplot(aes(x=run, y=value, fill=Key)) +
  geom_bar(stat="identity", position = "stack")+
  labs(x = "Sample", y="No. of reads predicted for all virus classes")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(angle=0,hjust=1),
        legend.position = "bottom"
        )+
  coord_flip()

p3 <- pcombdf %>% 
  select(run,class,perc_predicted_reads_95, classification,cutoff) %>%
  ggplot() +
  geom_point(aes(class,run, size=cutoff), color="black", alpha=0.3) +
  geom_point(aes(class,run, size=perc_predicted_reads_95,color=classification), alpha=1) +
  scale_color_manual(name="Classification", values=c("dodgerblue", "black"), breaks=c("Positive","Negative"))+
  labs(y = "Sample", x="Virus class")+
  guides(size=guide_legend(title="PR (%)", size=20))+
  #scale_x_discrete(labels= c("CoV","IV","hMPV", "hRV","SCoV2"))+
  theme_bw()+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme(
    axis.text.y=element_text(size=15),
    axis.text.x=element_text(angle=0,hjust=1,size=30),
    axis.title=element_blank(),
    legend.position = "bottom",
    legend.text=element_text(size=20),
    legend.title=element_text(size=20)
  )

p4 <-  pcombdf %>%
  select(run,classification,BWA_only,PACIFIC_only,overlap) %>%
  gather(overlap,BWA_only,PACIFIC_only, key="Key", value="value") %>%
  mutate(Key=factor(Key, levels=c("BWA_only", "PACIFIC_only", "overlap"))) %>%
  drop_na() %>% 
  ggplot(aes(x=run, y=value, fill=Key)) +
  geom_bar(stat="identity", position = "stack")+
  labs(x = "Sample", y="No. of reads predicted for all virus classes")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(angle=0,hjust=1,size=30),
        legend.position = "bottom",
  )+
  coord_flip()

#Create plots
plot_grid(p1, p2, rel_widths = c(2/3, 1/3))
plot_grid(p3, p4, rel_widths = c(2/3, 1/3))

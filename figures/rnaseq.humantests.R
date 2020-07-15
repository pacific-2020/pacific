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
#Cutoffs were determined as above 99% of the data in FPR tests using a beta distribution
co <- tibble(class = unique(dat$class), 
             cutoff = c(0.000213116230298516,0.0404740377229317,0.000806814626438465,0.000153651861475909,0.0418425381539617,NA,NA,NA)) %>% 
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
  mutate(run=filename,filename=paste(filename,"neg",sep="."))

pos <- semi_join(df,posmd,by=c("filename"="Run")) %>%
  mutate(run=filename,filename=paste(filename,"pos",sep="."))

#Join negative and positive tables
negpos <- rbind(neg,pos)
negpos$class <- gsub("SARS-CoV-2", "Sars_cov_2", negpos$class)

#Create bubble plot of classifications for positive and negative controls
negpos <- negpos %>% 
  select(filename,class,perc_predicted_reads_95, predicted_reads_95, classification,cutoff) %>%
  filter(!class %in% c("Human", "Discarded", "rc_discarded")) %>%
  arrange(desc(classification,filename)) %>%
  mutate(filename = fct_rev(as_factor(filename)))

#####
#Load read ids for further analysis
np <- read_tsv("pacific/figures/negpos.pacificfareadids.txt", col_names=F)
colnames(np) <-  c("run", "readid", "class")

ns <- read_tsv("pacific/figures/negposctrl.samreadids.unique.txt", col_names=F)
colnames(ns) <-  c("run", "readid", "class")

n20 <- read_tsv("pacific/figures/negposctrl.mapq20.samreadids.unique.txt", col_names=F)
colnames(n20) <-  c("run", "readid", "class")

nk <- read_tsv("pacific/figures/negpos_kraken_samples_out.txt", col_names=F)
nk <- nk %>% 
  select(X1,X2,X3)
colnames(nk) <-  c("run", "readid", "class")

#Text processing
np$run <- gsub("pacificoutput_", "", np$run)
ns$run <- gsub("\\..*", "", ns$run)
n20$run <- gsub("\\..*", "", n20$run)
ns$class <- gsub("Cornidovirineae", "Coronaviridae", ns$class)
n20$class <- gsub("Cornidovirineae", "Coronaviridae", n20$class)
nk$class <- gsub("sars_cov_2", "Sars_cov_2", nk$class)
nk$class <- gsub("Rhinovirus", "rhinovirus", nk$class)

np <- np %>% 
  mutate(id = paste(run,readid,class, sep=":"), pacific=class) %>%
  select(run,id,pacific)

ns <- ns %>% 
  mutate(id = paste(run,readid,class, sep=":"), bwa=class) %>%
  select(run,id,bwa)

n20 <- n20 %>% 
  mutate(id = paste(run,readid,class, sep=":"), bwa20=class) %>%
  select(run,id,bwa20)

nk <- nk %>% 
  mutate(id = paste(run,readid,class, sep=":"), kraken=class) %>%
  select(run,id,kraken)

nj <- full_join(np,ns, by=c("id","run"))
nj <- full_join(nj,n20, by=c("id","run"))
nj <- full_join(nj,nk, by=c("id","run"))

nj <- unique(nj) 

#write_tsv(nj, "/Users/Renzo/pbk.negpos.readidoverlap2.txt")

#Create summary numbers
pac <- nj %>%
  select(run, pacific) %>%
  group_by(run,pacific) %>%
  summarise(pac=n()) %>%
  rename("class"="pacific","pacific"="pac")

bw <- nj %>%
  select(run, bwa) %>%
  group_by(run,bwa) %>%
  summarise(b=n()) %>%
  rename("class"="bwa","bwa"="b")

bw20 <- nj %>%
  select(run, bwa20) %>%
  group_by(run,bwa20) %>%
  summarise(b20=n()) %>%
  rename("class"="bwa20","bwa20"="b20")

kra <- nj %>%
  select(run, kraken) %>%
  group_by(run,kraken) %>%
  summarise(kra=n()) %>%
  rename("class"="kraken","kraken"="kra")

tmp <- full_join(pac,bw, by=c("run","class"))
pbksum <- full_join(tmp,bw20, by=c("run","class"))
pbksum <- full_join(pbksum,kra, by=c("run","class"))

forj <- pbksum %>%
  select(run) %>%
  expand(class=c("Coronaviridae", "rhinovirus", "Sars_cov_2", "Metapneumovirus","Influenza"))

pbksummary <- left_join(forj,pbksum) %>%
  replace(is.na(.), 0)

pmd <- tibble(unique(pos$run),unique(pos$filename))
colnames(pmd) <- c("run","filename")
nmd <- tibble(unique(neg$run),unique(neg$filename))
colnames(nmd) <- c("run","filename")
pnmd <- rbind(pmd,nmd)

pbksummary <- left_join(pbksummary,pnmd)
pbksummary$filename <- factor(pbksummary$filename, levels=rev(unique(negpos$filename)),
)

pbksummary$class <- gsub("rhinovirus", "Rhinovirus", pbksummary$class)

#Read total read numbers
tot <- read_tsv("pacific/figures/negpossamples.totalreads.txt", col_names=F)
colnames(tot) <- c("run", "totalreads")
pbksummary <- left_join(pbksummary,tot)
pbksummary <- pbksummary %>%
  mutate(pacific_prop=(pacific/totalreads)*100,bwa_prop=(bwa/totalreads)*100,bwa20_prop=(bwa20/totalreads)*100,kraken_prop=(kraken/totalreads)*100, 
         cutoff = c(0.0404740377229317,0.000806814626438465,0.000153651861475909,0.0418425381539617,0.000213116230298516),
         pac_classification=ifelse(pacific_prop>=cutoff,"Positive","Negative"),
         bwa_classification=ifelse(bwa_prop>=cutoff,"Positive","Negative"),
         bwa20_classification=ifelse(bwa20_prop>=cutoff,"Positive","Negative"),
         krak_classification=ifelse(kraken_prop>=cutoff,"Positive","Negative"))

negpos.pbksummary <- pbksummary

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
kraken <- read_tsv("pacific/figures/gala.kraken.150.parsed.3col.txt", col_names=F)

colnames(fa) <-  c("run", "readid", "class")
colnames(sam) <-  c("run", "readid", "class")
colnames(kraken) <-  c("run", "readid", "class")

n20 <- read_tsv("pacific/figures/wastudy.mapq20.samreadids.unique.txt", col_names=F)
colnames(n20) <-  c("run", "readid", "class")
n20$run <- gsub("\\..*", "", n20$run)
n20$class <- gsub("Cornidovirineae", "Coronaviridae", n20$class)
kraken$class <- gsub("Rhinovirus", "rhinovirus", kraken$class)

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

krak <- kraken %>% 
  group_by(run,class) %>% 
  summarise(kraken=n())

bwa$class <- gsub("rhinovirus", "Rhinovirus", bwa$class )
krak$class <- gsub("rhinovirus", "Rhinovirus", krak$class )

#Prepare dataframes
s <- sam %>% 
  mutate(id = paste(run,readid,class, sep=":"), bwa=class) %>%
  select(run, id,bwa)

n20 <- n20 %>% 
  mutate(id = paste(run,readid,class, sep=":"), bwa20=class) %>%
  select(run,id,bwa20)

f <- fa %>% 
  mutate(id = paste(run,readid,class, sep=":"), pacific=class) %>%
  select(run,id,pacific)

k <- kraken %>% 
  mutate(id = paste(run,readid,class, sep=":"), kraken=class) %>%
  select(run,id,kraken)

j <- full_join(f,s, by=c("id","run"))
j <- full_join(j,n20, by=c("id","run"))
j <- full_join(j,k, by=c("id","run"))

j <- unique(j) 

#write_tsv(j, "pacific/metadata/pbk.wastudy.readidoverlap.txt")

#Create summary numbers
pac <- j %>%
  select(run, pacific) %>%
  group_by(run,pacific) %>%
  summarise(pac=n()) %>%
  rename("class"="pacific","pacific"="pac")

bw <- j %>%
  select(run, bwa) %>%
  group_by(run,bwa) %>%
  summarise(b=n()) %>%
  rename("class"="bwa","bwa"="b")

bw20 <- j %>%
  select(run, bwa20) %>%
  group_by(run,bwa20) %>%
  summarise(b20=n()) %>%
  rename("class"="bwa20","bwa20"="b20")

kra <- j %>%
  select(run, kraken) %>%
  group_by(run,kraken) %>%
  summarise(kra=n()) %>%
  rename("class"="kraken","kraken"="kra")

tmp <- full_join(pac,bw, by=c("run","class"))
pbksum <- full_join(tmp,bw20, by=c("run","class"))
pbksum <- full_join(pbksum,kra, by=c("run","class"))

forj <- pbksum %>%
  select(run) %>%
  expand(class=c("Coronaviridae", "rhinovirus", "Sars_cov_2", "Metapneumovirus","Influenza"))

pbksummary <- left_join(forj,pbksum) %>%
  replace(is.na(.), 0)

pbksummary$class <- gsub("rhinovirus", "Rhinovirus", pbksummary$class)

#Read total read numbers
tot <- read_tsv("pacific/figures/gala.totalreads.txt", col_names=F)
colnames(tot) <- c("run", "totalreads")
pbksummary <- left_join(pbksummary,tot)
pbksummary <- pbksummary %>%
  mutate(pacific_prop=(pacific/totalreads)*100,bwa_prop=(bwa/totalreads)*100,bwa20_prop=(bwa20/totalreads)*100,kraken_prop=(kraken/totalreads)*100, 
         cutoff = c(0.0404740377229317,0.000806814626438465,0.000153651861475909,0.0418425381539617,0.000213116230298516),
         pac_classification=ifelse(pacific_prop>=cutoff,"Positive","Negative"),
         bwa_classification=ifelse(bwa_prop>=cutoff,"Positive","Negative"),
         bwa20_classification=ifelse(bwa20_prop>=cutoff,"Positive","Negative"),
         krak_classification=ifelse(kraken_prop>=cutoff,"Positive","Negative")) %>%
  arrange(desc(pacific_prop,pac_classification, bwa_classification, kraken_classification)) %>%
  mutate(overallclass=ifelse(pac_classification =="Positive" | 
                               bwa_classification =="Positive" |
                               bwa20_classification =="Positive" |
                               krak_classification =="Positive", 
                             "VPositive", "VNegative")) %>%
  arrange(desc(overallclass, pac_classification,run, bwa_classification, kraken_classification))


pbksummary.all <- rbind(negpos.pbksummary,pbksummary)
pbksummary.all$overallclass <- factor(pbksummary.all$overallclass, levels = c("VPositive","VNegative"))
pbksummary.all$run <- factor(pbksummary.all$run, levels = rev(unique(pbksummary.all$run)))

#Draw plot for all positive samples
pacp <- pbksummary.all  %>%
  ggplot() +
  geom_point(aes(x=class,y=run, size=cutoff), color="black", alpha=0.3) +
  geom_point(aes(x=class,y=run, size=pacific_prop, color=pac_classification), alpha=1) +
  scale_color_manual(name="Classification", values=c("dodgerblue", "black"), breaks=c("Positive","Negative"))+
  labs(y = "Sample", x="Virus class")+
  guides(size=guide_legend(title="Read proportions (%)"))+
  theme_bw()+
  scale_size_continuous(breaks = c(0,0.2,0.4,0.6,0.8),
                       limits=c(0,1))+
  theme(
    axis.text.x=element_text(angle=30,hjust=1),
    axis.title=element_blank(),
    #text=element_text(size=16),
    #legend.title = element_blank(),
    legend.position = "bottom",
  )+
  ggtitle("PACIFIC")

bwap <- pbksummary.all %>%
  ggplot() +
  geom_point(aes(x=class,y=run, size=cutoff), color="black", alpha=0.3) +
  geom_point(aes(x=class,y=run, size=bwa_prop, color=bwa_classification), alpha=1) +
  scale_color_manual(name="Classification", values=c("dodgerblue", "black"), breaks=c("Positive","Negative"))+
  labs(y = "Sample", x="Virus class")+
  guides(size=guide_legend(title="PR (%)"))+
  theme_bw()+
  scale_size_continuous(breaks = c(0,0.2,0.4,0.6,0.8),
                        limits=c(0,1))+
  theme(
    axis.text.x=element_text(angle=30,hjust=1),
    axis.text.y=element_blank(),
    axis.title=element_blank(),
    #text=element_text(size=16),
    legend.title = element_blank(),
    legend.position = "bottom",
  )+
  ggtitle("BWA")

bwa20p <- pbksummary.all %>%
  ggplot() +
  geom_point(aes(x=class,y=run, size=cutoff), color="black", alpha=0.3) +
  geom_point(aes(x=class,y=run, size=bwa20_prop, color=bwa20_classification), alpha=1) +
  scale_color_manual(name="Classification", values=c("dodgerblue", "black"), breaks=c("Positive","Negative"))+
  labs(y = "Sample", x="Virus class")+
  guides(size=guide_legend(title="PR (%)"))+
  theme_bw()+
  scale_size_continuous(breaks = c(0,0.2,0.4,0.6,0.8),
                        limits=c(0,1))+
  theme(
    axis.text.x=element_text(angle=30,hjust=1),
    axis.text.y=element_blank(),
    axis.title=element_blank(),
    #text=element_text(size=16),
    legend.title = element_blank(),
    legend.position = "bottom",
  )+
  ggtitle("BWA >MAPQ 20")

krakp <- pbksummary.all %>%
  ggplot() +
  geom_point(aes(x=class,y=run, size=cutoff), fill="black", alpha=0.3) +
  geom_point(aes(x=class,y=run, size=kraken_prop, color=krak_classification), alpha=1) +
  scale_color_manual(name="Classification", values=c("dodgerblue", "black"), breaks=c("Positive","Negative"))+
  labs(y = "Sample", x="Virus class")+
  guides(size=guide_legend(title="Number of reads"))+
  theme_bw()+
  scale_size_continuous(breaks = c(0,0.2,0.4,0.6,0.8),
                        limits=c(0,1))+
  theme(
    axis.text.x=element_text(angle=30,hjust=1),
    axis.text.y=element_blank(),
    #text=element_text(size=16),
    axis.title=element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
  )+
  ggtitle("Kraken2")

plot_grid(pacp,bwap, bwa20p, krakp, rel_widths = c(2.5/6,1.75/6, 1.75/6, 1.75/6), ncol=4)

#Draw plot for positive samples only
positive <- negpos[negpos$classification %in% "Positive",]
positive$filename <- gsub(".pos", "", positive$filename)
positive <- negpos.pbksummary %>% 
  filter(run %in% positive$filename)
positivewa <- pbksummary[pbksummary$overallclass %in% "VPositive",]
positivewa <- pbksummary %>% 
  filter(run %in% positivewa$run)
positiveall <- rbind(positive,positivewa)
positiveall$run <- factor(positiveall$run, levels = rev(unique(positiveall$run)))

#Draw plot for all positive samples
pacp <- positiveall  %>%
  ggplot() +
  geom_point(aes(x=class,y=run, size=cutoff), color="black", alpha=0.3) +
  geom_point(aes(x=class,y=run, size=pacific_prop, color=pac_classification), alpha=1) +
  scale_color_manual(name="Classification", values=c("dodgerblue", "black"), breaks=c("Positive","Negative"))+
  labs(y = "Sample", x="Virus class")+
  guides(size=guide_legend(title="Read proportions (%)"))+
  theme_bw()+
  scale_size_continuous(breaks = c(0,0.2,0.4,0.6,0.8),
                        limits=c(0,1))+
  theme(
    axis.text.x=element_text(angle=30,hjust=1),
    axis.title=element_blank(),
    #text=element_text(size=16),
    #legend.title = element_blank(),
    legend.position = "bottom",
  )+
  ggtitle("PACIFIC")

bwap <- positiveall %>%
  ggplot() +
  geom_point(aes(x=class,y=run, size=cutoff), color="black", alpha=0.3) +
  geom_point(aes(x=class,y=run, size=bwa_prop, color=bwa_classification), alpha=1) +
  scale_color_manual(name="Classification", values=c("dodgerblue", "black"), breaks=c("Positive","Negative"))+
  labs(y = "Sample", x="Virus class")+
  guides(size=guide_legend(title="PR (%)"))+
  theme_bw()+
  scale_size_continuous(breaks = c(0,0.2,0.4,0.6,0.8),
                        limits=c(0,1))+
  theme(
    axis.text.x=element_text(angle=30,hjust=1),
    axis.text.y=element_blank(),
    axis.title=element_blank(),
    #text=element_text(size=16),
    legend.title = element_blank(),
    legend.position = "bottom",
  )+
  ggtitle("BWA")

krakp <- positiveall %>%
  ggplot() +
  geom_point(aes(x=class,y=run, size=cutoff), fill="black", alpha=0.3) +
  geom_point(aes(x=class,y=run, size=kraken_prop, color=krak_classification), alpha=1) +
  scale_color_manual(name="Classification", values=c("dodgerblue", "black"), breaks=c("Positive","Negative"))+
  labs(y = "Sample", x="Virus class")+
  guides(size=guide_legend(title="Number of reads"))+
  theme_bw()+
  scale_size_continuous(breaks = c(0,0.2,0.4,0.6,0.8),
                        limits=c(0,1))+
  theme(
    axis.text.x=element_text(angle=30,hjust=1),
    axis.text.y=element_blank(),
    #text=element_text(size=16),
    axis.title=element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
  )+
  ggtitle("Kraken2")

plot_grid(pacp,bwap, krakp, rel_widths = c(2.75/6,1.75/6, 1.75/6), ncol=3)

#Extract 10 negative samples for final experiment
negative <- negpos[negpos$classification %in% "Negative",]
negative$filename <- gsub(".neg", "", negative$filename)
negative <- negpos.pbksummary %>% 
  filter(run %in% negative$filename) %>%
  arrange(desc(bwa_classification))
negative$run <- factor(negative$run, levels = rev(unique(negative$run)))

#Draw plot for negative samples
pacp <- negative  %>%
  ggplot() +
  geom_point(aes(x=class,y=run, size=cutoff), color="black", alpha=0.3) +
  geom_point(aes(x=class,y=run, size=pacific_prop, color=pac_classification), alpha=1) +
  scale_color_manual(name="Classification", values=c("dodgerblue", "black"), breaks=c("Positive","Negative"))+
  labs(y = "Sample", x="Virus class")+
  guides(size=guide_legend(title=""))+
  theme_bw()+
  scale_size_continuous(breaks = c(0,0.1,0.2,0.3,0.4),
                        limits=c(0,0.4))+
  theme(
    axis.text.x=element_text(angle=30,hjust=1),
    axis.title=element_blank(),
    #text=element_text(size=16),
    #legend.title = element_blank(),
    legend.position = "bottom",
  )+
  ggtitle("PACIFIC")

bwap <- negative %>%
  ggplot() +
  geom_point(aes(x=class,y=run, size=cutoff), color="black", alpha=0.3) +
  geom_point(aes(x=class,y=run, size=bwa_prop, color=bwa_classification), alpha=1) +
  scale_color_manual(name="Classification", values=c("dodgerblue", "black"), breaks=c("Positive","Negative"))+
  labs(y = "Sample", x="Virus class")+
  guides(size=guide_legend(title="PR (%)"))+
  theme_bw()+
  scale_size_continuous(breaks = c(0,0.1,0.2,0.3,0.4),
                        limits=c(0,0.4))+
  theme(
    axis.text.x=element_text(angle=30,hjust=1),
    axis.text.y=element_blank(),
    axis.title=element_blank(),
    #text=element_text(size=16),
    legend.title = element_blank(),
    legend.position = "bottom",
  )+
  ggtitle("BWA")

krakp <- negative %>%
  ggplot() +
  geom_point(aes(x=class,y=run, size=cutoff), fill="black", alpha=0.3) +
  geom_point(aes(x=class,y=run, size=kraken_prop, color=krak_classification), alpha=1) +
  scale_color_manual(name="Classification", values=c("dodgerblue", "black"), breaks=c("Positive","Negative"))+
  labs(y = "Sample", x="Virus class")+
  guides(size=guide_legend(title="Number of reads"))+
  theme_bw()+
  scale_size_continuous(breaks = c(0,0.1,0.2,0.3,0.4),
                        limits=c(0,0.4))+
  theme(
    axis.text.x=element_text(angle=30,hjust=1),
    axis.text.y=element_blank(),
    #text=element_text(size=16),
    axis.title=element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
  )+
  ggtitle("Kraken2")

plot_grid(pacp,bwap, krakp, rel_widths = c(2.75/6,2/6, 2/6), ncol=3)

#Create metadata table
plotdata <- data.frame()
plotdata <- bind_rows(plotdata, select(pbksummary.all, run,class,cutoff,nreads=pacific,rproportion=pacific_prop,classification=pac_classification, overallclass) %>% mutate(classifier="pacific"))
plotdata <- bind_rows(plotdata, select(pbksummary.all, run,class,cutoff,nreads=bwa,rproportion=bwa_prop,classification=bwa_classification, overallclass) %>% mutate(classifier="bwa"))
plotdata <- bind_rows(plotdata, select(pbksummary.all, run,class,cutoff, nreads=kraken,rproportion=kraken_prop,classification=krak_classification, overallclass) %>% mutate(classifier="kraken"))

options(scipen = 999)
plotdata <- plotdata %>%
  select(run,classifier,class,nreads,rproportion) %>%
  arrange(as.character(run))

pmd <- pmd %>%
  mutate(study="SARS-CoV-2+_technical_replicate") %>%
  arrange(as.character(run))
pmd[1,]$study <- "SARS-CoV-2+_human"

nmd <- nmd %>%
  mutate(study="Negative")%>%
  arrange(as.character(run))

gala <- tibble(run=unique(pbksummary$run))
gala <- gala %>% 
  mutate(filename=run,study="GALA") %>%
  arrange(as.character(run))

allmd <- rbind(pmd,nmd) %>%
  rbind(.,gala) %>%
  select(run,study)

finalmd <- left_join(plotdata,allmd, by="run") %>%
  arrange(as.character(run))

write_tsv(finalmd, "pacific/metadata/rnaseq.humantests.txt")

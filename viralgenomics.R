library(tidyverse)
library(rentrez)
library(xml2)
library(XML)
library(rbenchmark)
library(compiler)
library(RSQLite)
library(Biobase)
library(RCurl)
library(ggthemes)
library(ggrepel)

unescape_html <- function(str){
  xt <- xml_text(read_html(paste0("<x>", str, "</x>")))
  return(paste("<meta>", xt, "</meta>"))
}

simplifycolumn <- function(x) {
  return(sapply(x, paste, collapse =":"))
}

getSummaryInfo <- function(iteration, webhistory, maxrecords) {
  print(paste("processing", iteration, "iteration"))
  sumInfo <- entrez_summary("assembly",web_history = webhistory, retstart = iteration*maxrecords+1, retmax=maxrecords)
  sumtable <- data.frame(matrix(extract_from_esummary(sumInfo, names(sumInfo[[1]])), nrow = length(sumInfo), byrow = T))
  colnames(sumtable) <- names(sumInfo[[1]])
  return(sumtable)
}

processmeta <- function(x) {
  #a <- as_list(read_xml(unescape_html(x)))
  a <- xmlToList(xmlParse(paste("<doc>",unescape_html(x),"</doc>"), asText = T))
  a <- a$meta[["Stats"]]
  a <- as.data.frame(t(matrix(unlist(a), nrow = 3)))
  row.names(a) <- paste(a$V2,a$V3,sep=".")
  a <- select(a,V1)
  return(a)
}

processmeta <- cmpfun(processmeta)
getSummaryInfo <- cmpfun(getSummaryInfo)
simplifycolumn <- cmpfun(simplifycolumn)
unescape_html <- cmpfun(unescape_html)

##get a list of viral genome assembly from NCBI
sinfo <- NULL
viralgenomes <- entrez_search("assembly", 'txid10239[Organism:exp] & latest[filter] & "complete genome"[filter]', use_history = T)
maxrecords <- 500
iterations <- as.integer(viralgenomes$count/maxrecords)
start <- nrow(sinfo)/maxrecords

for (i in 0:iterations) {
  tmp <- getSummaryInfo(iteration = i, webhistory = viralgenomes$web_history, maxrecords=maxrecords)
  sinfo <- bind_rows(sinfo, tmp)
}

y <- sinfo
##property list is not ordered so just paste it together
y$propertylist <- sapply(y$propertylist, paste, collapse =":")

##process meta column to extract stats information out of it
x <- do.call("cbind",lapply(y$meta, processmeta))
x <- t(x)
x <- type.convert(x)
y <- cbind(y, x)
y$meta <- unlist(lapply(y$meta, unescape_html)) ##need to sanitise this. runs twice

##process multi item columns
y <- y %>% unnest_wider(col = synonym, names_sep = ".", names_repair = "universal")
y <- y %>% unnest_wider(col = anomalouslist, names_sep = ".", names_repair = "universal")
y <- y %>% unnest_wider(col = biosource, names_sep = ".", names_repair = "universal")
y <- y %>% unnest_wider(col = gb_bioprojects, names_sep = ".", names_repair = "universal")
y <- y %>% unnest_wider(col = rs_bioprojects, names_sep = ".", names_repair = "universal")
y <- y %>% unnest_wider(col = biosource.infraspecieslist, names_sep = ".", names_repair = "universal")
##unlist column items lists
y <- data.frame(apply(y,2, simplifycolumn))
##remove sortorder column
y <- select(y, -sortorder, -biosource.infraspecieslist.)
##do type convert for columns
y <- type.convert(y)
##date columns are not formatted correctly using type convert so take care of them
y <- mutate_at(y, colnames(y)[grep("date", colnames(y))], as.Date)


##when dates are written in sqlite table, they are written as real numbers
##to convert them back to legible date format
##as.Date(18354.0, "1970-01-01") ##origin must be specified
##Sys.Date() as todays date
##taxid should be used instead of speciestaxid. 
##(1) There are more unique entries in taxid vs speciestaxid
##(2) TaxID seems to refer to a strain level identification vs speciestaxid as name suggests
##There can be more than one assembly per taxid, probably because of different strains

##get Lineage information
download.file("https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz", "/media/dstore/covid19/ncbitaxdmp/new_taxdump.tar.gz")
setwd("/media/dstore/covid19/ncbitaxdmp/")
system("tar -xvzf new_taxdump.tar.gz")
##lineage information
system('sed "s/\t//g" rankedlineage.dmp | sed "s/\'//g" >rankedlineage.dmp.tmp')
rankedlineage <- read.table("rankedlineage.dmp.tmp", sep="|", quote = "", fill = T)
colnames(rankedlineage) <- c("taxid", "tax_name", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom", "empty")
##host information
system('sed "s/\t//g" host.dmp >host.dmp.tmp')
host <- read.table("host.dmp.tmp", sep="|", quote ="")
colnames(host) <- c("taxid", "potentialhosts")
##merged lineage information

z <- merge(y, rankedlineage, by.x = "taxid", by.y = "taxid", all.x = T)
z <- merge(z, host, by.x = "taxid", by.y = "taxid", all.x = T)
z <- droplevels(z)
z <- select(z, -empty)




##consolidate host information
z <- z %>% mutate(potentialhosts = case_when( .$potentialhosts == "vertebrates,human" ~ "human,vertebrates", is.na(.$potentialhosts) ~ "Unknown", .$potentialhosts == "vertebrates,invertebrates,human" ~ "human,invertebrates,vertebrates", .$potentialhosts == "invertebrates,vertebrates" ~ "vertebrates,invertebrates", TRUE ~ as.character(.$potentialhosts))) %>% mutate(potentialhosts = as.factor(potentialhosts))

##identify representative genome for a given tax id
spectab <- data.frame(table(z$taxid))
colnames(spectab) <- c("taxid", "genomes")
spectab$taxid <- as.integer(as.character(spectab$taxid))

##create sqlite database for later use
db = dbConnect(SQLite(), dbname="/media/dstore/covid19/ncbi.virus.20200409.sqlite")
dbWriteTable(db, "assembly", z)



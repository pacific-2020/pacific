library(tidyverse)
library(rentrez)
library(xml2)
library(XML)
library(rbenchmark)
library(compiler)
library(RSQLite)
library(Biobase)

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

##create sqlite database for later use
db = dbConnect(SQLite(), dbname="/media/dstore/covid19/ncbi.virus2.sqlite")
dbWriteTable(db, "assembly", y)

##when dates are written in sqlite table, they are written as real numbers
##to convert them back to legible date format
##as.Date(18354.0, "1970-01-01") ##origin must be specified
##Sys.Date() as todays date
##taxid should be used instead of speciestaxid. 
##(1) There are more unique entries in taxid vs speciestaxid
##(2) TaxID seems to refer to a strain level identification vs speciestaxid as name suggests
##There can be more than one assembly per taxid, probably because of different strains


#viraltaxonomy <- entrez_search("taxonomy", 'txid10239[Organism:exp]', use_history = T)

virustaxid <- unique(y$taxid)
getLineage <- function(start, maxrecords) {
  end <- start + maxrecords - 1
  end <- ifelse(end>length(virustaxid), length(virustaxid), end)
  taxml <- entrez_fetch("taxonomy",id = virustaxid[start:end], rettype = "xml")
  taxlineage <- xmlToList(xmlParse(taxml, asText = T))
  l <- data.frame(taxid = unlist(subListExtract(taxlineage, 'TaxId')), lineage = unlist(subListExtract(taxlineage, 'Lineage')))
  return(l)
}
maxrecords <- 200
lineageinfo <- do.call("rbind", lapply(seq(1,length(virustaxid),maxrecords), getLineage, maxrecords))

#y <- merge(x, lineageinfo, by.x = "taxid", by.y = "taxid", all.x = T)
#seqtable <- read.table("/media/dstore/covid19/seqid.vs.genomeid.txt", header = F, sep = "\t")
#z <- merge(seqtable, y, by.x = "V1", by.y = "synonym_gbacc", all.x = T)
#write.table(z, "/media/dstore/covid19/seqid.vs.genomeinfo.txt", sep = "\t", col.names = T, row.names = F, quote = F)



#entrez_fetch("assembly", web_history = viralgenomes$web_history, retstart = 1, retmax = 20, rettype = "xml", retmode = "full")

#getSummaryInfo <- function(iteration, webhistory, maxrecords) {
#  print(paste("processing", iteration, "iteration"))
#  x <- read_xml(entrez_fetch("assembly", web_history = webhistory, retstart = iteration*maxrecords+1, retmax = maxrecords, rettype = "docsum", retmode = "xml"))
#  y <- bind_rows(as_list(x))
#  return(y)
#}

#viralgenomes <- entrez_search("assembly", 'txid10239[Organism:exp] & latest[filter] & "complete genome"[filter]', use_history = T)
#maxrecords <- 200
#iterations <- as.integer(viralgenomes$count/maxrecords)
#sinfo <- do.call("bind_rows", lapply(c(0:iterations), getSummaryInfo, webhistory = viralgenomes$web_history, maxrecords=maxrecords))


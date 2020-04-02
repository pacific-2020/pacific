library(rentrez)
library(XML)
library(xml2)
library(tidyverse)
library(dplyr)
library(rbenchmark)
library(compiler)
library(stringr)
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
viralgenomes <- entrez_search("assembly", 'txid10239[Organism:exp] & latest[filter] & "complete genome"[filter]', use_history = T)
maxrecords <- 500
iterations <- as.integer(viralgenomes$count/maxrecords)

sinfo <- do.call("rbind", lapply(c(0:iterations), getSummaryInfo, webhistory = viralgenomes$web_history, maxrecords=maxrecords))
csinfo <- data.frame(apply(sinfo, 2, simplifycolumn))

numcol <- c("uid", "rsuid", "gbuid", "chainid", "taxid", "speciestaxid", "coverage", "primary", "contign50", "scaffoldn50")
factorcol <- c("assemblyaccession", "lastmajorreleaseaccession", "latestaccession", "assemblyname", "ucscname", "ensemblname", "speciesname", "assemblytype", "assemblyclass", "assemblystatus","biosampleaccn","biosampleid", "releaselevel", "releasetype", "refseq_category", "exclfromrefseq")
boolcol <- c("wgs", "gb_projects", "rs_projects", "assemblydescription")
charcol <- c("organism", "gb_bioprojects", "rs_bioprojects", "biosource", "partialgenomerepresentation", "submitterorganization", "anomalouslist", "propertylist", "fromtype", "synonym", "ftppath_genbank", "ftppath_refseq", "ftppath_assembly_rpt", "ftppath_stats_rpt", "ftppath_regions_rpt", "sortorder", "meta")
datecol <- c("asmreleasedate_genbank","asmreleasedate_refseq","seqreleasedate", "asmupdatedate", "submissiondate", "lastupdatedate")

csinfo[,numcol] <- apply(csinfo[,numcol], 2, function(x) as.numeric(as.character(x)))
csinfo[,charcol] <- lapply(csinfo[,charcol], as.character)
csinfo[,datecol] <- lapply(csinfo[,datecol], as.Date)


##partialgenomerepresentation values are char string of false
##anomalouslist is a list of zero or two elements
##propertylist can have anywhere from 6 to 11 items. Perhaps later sort it to unify it.
#cols = c(1, 3, 4, 5);    
#df[,cols] = apply(df[,cols], 2, function(x) as.numeric(as.character(x)));

##gb_bioprojects is a list of two elements stitched by ':', BioprojectAccn and BioprojectId
##rs_bioprojects is a list of two elements stitched by ':', BioprojectAccn and BioprojectId
csinfo <- mutate(csinfo, gb_bioprojectaccn = sapply(strsplit(csinfo$gb_bioprojects,":"), `[`, 1))
csinfo <- mutate(csinfo, gb_bioprojectid = sapply(strsplit(csinfo$gb_bioprojects,":"), `[`, 2))
csinfo$gb_bioprojectid <- as.numeric(as.character(csinfo$gb_bioprojectid))
csinfo <- mutate(csinfo, rs_bioprojectaccn = sapply(strsplit(csinfo$rs_bioprojects,":"), `[`, 1))
csinfo <- mutate(csinfo, rs_bioprojectid = sapply(strsplit(csinfo$rs_bioprojects,":"), `[`, 2))
csinfo$rs_bioprojectid <- as.numeric(as.character(csinfo$rs_bioprojectid))
##synonym is a combination of three values, GenBank Acc, RefSeqAcc, and Identical. If refseq and genbank are identical, it will say identical in the third column
csinfo <- mutate(csinfo, synonym_gbacc = sapply(strsplit(csinfo$synonym,":"), `[`, 1))
csinfo <- mutate(csinfo, synonym_rsacc = sapply(strsplit(csinfo$synonym,":"), `[`, 2))
csinfo <- mutate(csinfo, synonym_identical = sapply(strsplit(csinfo$synonym,":"), `[`, 3))
##biosource is a list of three elements stitched ':', InfraspeciesList, Sex and Isolate. InfraspeciesList is a list of two, Infraspecie sub_type and sub_value
m <- str_match(csinfo$biosource, "^list\\((.*?)\\):(.*?):(.*)")
m[,2] <- gsub("sub_type = \"", "", m[,2])
m[,2] <- gsub("\", sub_value = \""," = ",m[,2])
m[,2] <- gsub("\"$","",m[,2])
csinfo <- mutate(csinfo, biosource_strain = m[,2])
csinfo <- mutate(csinfo, biosource_sex = m[,3])
csinfo <- mutate(csinfo, biosource_isolate = m[,4])


##process meta column
##see https://github.com/hp2048/teaching/blob/master/assemblyInfo.R for details if required
x <- do.call("cbind",lapply(csinfo$meta, processmeta))
x <- t(x)
x <- type.convert(x)
csinfo <- cbind(csinfo, x)
csinfo <- select(csinfo, -meta, -sortorder, -gb_bioprojects, -rs_bioprojects, -synonym, -biosource)

write.table(csinfo, "/media/dstore/covid19/assemblyInfo.20200402.txt", sep = "\t", quote = F, row.names = F, col.names = T)




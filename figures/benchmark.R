library(Rmisc)
library(tidyverse)
library(grid)
library(gridExtra)
library(RColorBrewer)


calcmetrics <- function(a, datatype, givenlabel) {
  ##metrics calculated according to definitions at https://en.wikipedia.org/wiki/Precision_and_recall
  pdiscarded <- sum(a %>% filter(knownlabel == givenlabel & predictedlabel == "pDiscarded") %>% .$readcounts)
  cdiscarded <- sum(a %>% filter(knownlabel == givenlabel & predictedlabel == "rc_discarded") %>% .$readcounts)
  discarded  <- sum(a %>% filter(knownlabel == givenlabel & predictedlabel == "Discarded") %>% .$readcounts)
  totalreads <- sum(a %>% filter(knownlabel == givenlabel) %>% .$readcounts)
  #a <- a %>% filter(predictedlabel != "Discarded" & predictedlabel != "pDiscarded")
  a <- filter(a, predictedlabel != "Discarded")
  truepositive  <- sum( a %>% filter(predictedlabel == givenlabel & knownlabel == givenlabel) %>% .$readcounts)
  falsepositive <- sum( a %>% filter(predictedlabel == givenlabel & knownlabel != givenlabel) %>% .$readcounts)
  truenegative  <- sum( a %>% filter(predictedlabel != givenlabel & knownlabel != givenlabel) %>% .$readcounts)
  falsenegative <- sum( a %>% filter(predictedlabel != givenlabel & knownlabel == givenlabel) %>% .$readcounts)

  precision <- truepositive/(truepositive + falsepositive)
  recall    <- truepositive/(truepositive + falsenegative)
  accuracy  <- (truepositive + truenegative) / (truepositive + falsepositive + truenegative + falsenegative)
  tpr       <- truepositive / (truepositive + falsenegative)
  fpr       <- falsepositive / (falsepositive + truenegative)
  falsediscoveryrate <- falsepositive / (falsepositive + truepositive)
  tnr       <- truenegative / (truenegative + falsepositive)
  fnr       <- falsenegative / (falsenegative + truepositive)
  falseomissionrate <- falsenegative / (falsenegative + truenegative)
  f1 <- 2 * ((precision * tpr) / (precision + tpr))
  plr <- tpr/fpr
  nlr <- fnr/tnr
  dor <- plr/nlr
  prerec <- data.frame(dataid = did,
                       datatype = datatype,
                       class = givenlabel,
                       total = totalreads,
                       tp = truepositive,
                       fp = falsepositive,
                       tpr = tpr,
                       fpr = fpr,
                       falsediscoveryrate = falsediscoveryrate,
                       tn = truenegative,
                       fn = falsenegative,
                       tnr = tnr,
                       fnr = fnr,
                       falseomissionrate = falseomissionrate,
                       precision = precision,
                       recall = recall,
                       accuracy = accuracy,
                       balancedaccuracy = (tpr + tnr)/2,
                       f1 = f1,
                       plr = plr,
                       nlr = nlr,
                       dor = dor,
                       discarded = discarded,
                       pdiscarded = pdiscarded,
                       cdiscarded = cdiscarded
  )
  prerec
}

predictions <- read.table("benchmarkresults.txt", sep = "\t", header = T)
predictionmetrics <- data.frame()
for (did in paste("T",seq(1:100),sep="")){
  for (givenlabel in levels(predictions$knownlabel)) {
    for (pt in levels(predictions$ptype)) {
      a <- predictions %>% mutate_if(is.factor, as.character) %>% filter(dataid == did & ptype == pt)
      predictionmetrics <- bind_rows(predictionmetrics, calcmetrics(a = a, datatype = "all", givenlabel = givenlabel) %>% mutate(ptype = pt)) 
      b <- predictions %>% mutate_if(is.factor, as.character) %>% filter( dataid == did & readtype != "EQUAL" & ptype == pt)
      predictionmetrics <- bind_rows(predictionmetrics, calcmetrics(a = b, datatype = "mm", givenlabel = givenlabel) %>% mutate(ptype = pt))
      c <- predictions %>% mutate_if(is.factor, as.character) %>% filter( dataid == did & readtype == "EQUAL" & ptype == pt)
      predictionmetrics <- bind_rows(predictionmetrics, calcmetrics(a = c, datatype = "equal", givenlabel = givenlabel) %>% mutate(ptype = pt))
    }
  }
}

pm <- predictionmetrics
predictionmetrics <- predictionmetrics %>% mutate(class = case_when(class=="Sars_cov_2" ~ "SARS-CoV-2", TRUE ~ as.character(class))) %>% mutate(class=as.factor(class))


###Table 1: PACIFIC performance metrics for each class in 100 independent simulated datasets
metricstable <- data.frame()
for (m in c("fnr", "fpr", "precision", "recall", "accuracy")) {
  metricstable <- bind_rows(metricstable,
                            filter(predictionmetrics, datatype == "all" & ptype == "Double" & class != "Unrelated") %>% 
                              summarySE(measurevar = m, groupvars = c("class")) %>% 
                              dplyr::select(class,all_of(m),ci) %>% 
                              mutate_if(is.numeric, ~format(., nsmall = 1)) %>% 
                              mutate(test = m, op="(",cp=")",s=" ") %>% 
                              unite(val, m, s, op,ci,cp, sep = "")
  )
}

spread(metricstable, test, val)
write.table(spread(metricstable, test, val), file="performance.metrics.txt", sep = "\t", col.names = T, row.names = F, quote = F)

###Figure 2:
grid_arrange_shared_legend <-
  function(...,
           ncol = length(list(...)),
           nrow = 1,
           position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
  }

fprplot <- predictionmetrics %>% 
  dplyr::filter(datatype=="all" & ! class %in% c("Human", "Unrelated")) %>% 
  dplyr::mutate(ptype=factor(ptype, levels = c("Single", "Double"))) %>% 
  ggplot(aes(x=class,y=fpr,fill=ptype)) + 
  geom_boxplot(notch = T, lwd = 0.3) +
  labs(fill = "Prediction type") +
  scale_fill_discrete(labels = c("Single", "Double")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), text = element_text(size=16)) +
  ylab("False positive rate") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  xlab("Virus class") 

fnrplot <- predictionmetrics %>% 
  dplyr::filter(datatype=="all" & ! class %in% c("Human", "Unrelated")) %>% 
  dplyr::mutate(ptype=factor(ptype, levels = c("Single", "Double"))) %>% 
  ggplot(aes(x=class,y=fnr,fill=ptype)) + 
  geom_boxplot(notch = T, lwd = 0.3) +
  labs(fill = "Prediction type") +
  scale_fill_discrete(labels = c("Single", "Double")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), text = element_text(size=16)) +
  ylab("False negative rate") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  xlab("Virus class") 


grid_arrange_shared_legend(fprplot,fnrplot, nrow = 1)


###Figure 3:
### comparisons of FPR and FNR for read types
fpr3 <- predictionmetrics %>% 
  filter(ptype == "Double" & datatype != "all" & class != "Human" & class != "Unrelated") %>% 
  dplyr::select(dataid, class, datatype, fpr, fnr) %>% 
  gather("metric", "value", -dataid, -class, -datatype) %>% 
  mutate(metric = case_when(metric == "fpr" ~ "False positive rate", metric == "fnr" ~ "False negative rate", TRUE ~ as.character(metric))) %>% 
  mutate(metric = factor(metric, levels = c("False positive rate", "False negative rate"))) %>%
  mutate(datatype = case_when(datatype == "mm" ~ "Mismatch", datatype == "equal" ~ "Exact", TRUE ~ as.character(datatype)), datatype = as.factor(datatype)) %>%
  dplyr::filter(metric == "False positive rate") %>%
  ggplot(aes(x=class,y=value, fill = datatype)) + 
  geom_boxplot(notch = T, lwd = 0.3) + 
  theme_bw()+ 
  xlab("Virus class") + ylab("False positive rate") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  labs(fill = "Read types") + theme(axis.text.x = element_text(angle = 30, hjust = 1), text = element_text(size=16))

fnr3 <- predictionmetrics %>% 
  filter(ptype == "Double" & datatype != "all" & class != "Human" & class != "Unrelated") %>% 
  dplyr::select(dataid, class, datatype, fpr, fnr) %>% 
  gather("metric", "value", -dataid, -class, -datatype) %>% 
  mutate(metric = case_when(metric == "fpr" ~ "False positive rate", metric == "fnr" ~ "False negative rate", TRUE ~ as.character(metric))) %>% 
  mutate(metric = factor(metric, levels = c("False positive rate", "False negative rate"))) %>%
  mutate(datatype = case_when(datatype == "mm" ~ "Mismatch", datatype == "equal" ~ "Exact", TRUE ~ as.character(datatype)), datatype = as.factor(datatype)) %>%
  dplyr::filter(metric == "False negative rate") %>%
  ggplot(aes(x=class,y=value, fill = datatype)) + 
  geom_boxplot(notch = T, lwd = 0.3) + 
  theme_bw()+ 
  xlab("Virus class") + ylab("False negative rate") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  labs(fill = "Read types") + theme(axis.text.x = element_text(angle = 30, hjust = 1), text = element_text(size=16))

grid_arrange_shared_legend(fpr3,fnr3, nrow = 1)


###Figure 4:
### Sources of FP and FN
###false positive sources bar chart
#mycolors <- sample(brewer.pal(n=12,name="Paired"), size = 7, replace = F )
#print(mycolors)
p1 <- predictions %>% 
  dplyr::mutate(knownlabel = case_when(knownlabel=="Sars_cov_2" ~ "SARS-CoV-2", TRUE ~ as.character(knownlabel))) %>% 
  dplyr::mutate(knownlabel=as.factor(knownlabel)) %>%
  dplyr::mutate(predictedlabel = case_when(predictedlabel=="Sars_cov_2" ~ "SARS-CoV-2", TRUE ~ as.character(predictedlabel))) %>% 
  dplyr::mutate(predictedlabel=as.factor(predictedlabel)) %>%
  dplyr::filter(as.character(knownlabel) != as.character(predictedlabel) & ptype=="Double") %>% 
  dplyr::group_by(dataid,knownlabel,predictedlabel) %>% 
  dplyr::summarise(totalcounts = sum(readcounts)) %>% 
  dplyr::group_by(knownlabel, predictedlabel) %>% 
  dplyr::summarise(average = sum(totalcounts)/100) %>% 
  dplyr::filter(! predictedlabel %in% c("Discarded", "Human", "pDiscarded", "rc_discarded")) %>% 
  ggplot(aes(x=predictedlabel, y = average, fill = knownlabel)) + 
  geom_bar(stat="identity", position = "stack") + 
  ylab("Mean false positive reads") + 
  xlab("Predicted class") +
  theme_bw() + labs(fill = "True\nClass") +
  theme(text = element_text(size=12), axis.text.x = element_text(angle = 30, hjust = 1), legend.position="bottom") +
  scale_fill_manual(values = c("#33A02C", "#B2DF8A" ,"#A6CEE3" ,"#FDBF6F", "#E31A1C" ,"#FB9A99", "#1F78B4"))

p2 <- predictions %>% 
  dplyr::mutate(knownlabel = case_when(knownlabel=="Sars_cov_2" ~ "SARS-CoV-2", TRUE ~ as.character(knownlabel))) %>% 
  dplyr::mutate(knownlabel=as.factor(knownlabel)) %>%
  dplyr::mutate(predictedlabel = case_when(predictedlabel=="Sars_cov_2" ~ "SARS-CoV-2", TRUE ~ as.character(predictedlabel))) %>% 
  dplyr::mutate(predictedlabel=as.factor(predictedlabel)) %>%
  dplyr::filter(as.character(knownlabel) != as.character(predictedlabel) & ptype=="Double") %>% 
  dplyr::group_by(dataid,knownlabel,predictedlabel) %>% 
  dplyr::summarise(totalcounts = sum(readcounts)) %>% 
  dplyr::group_by(knownlabel, predictedlabel) %>% 
  dplyr::summarise(average = sum(totalcounts)/100) %>% 
  dplyr::filter(predictedlabel != "Discarded" & average > 0 & knownlabel != "Unrelated" & knownlabel != "Human") %>% 
  ggplot(aes(x=knownlabel, y = average, fill = predictedlabel)) + 
  geom_bar(stat="identity", position = "stack") + 
  ylab("Mean false negative reads") + 
  xlab("True class") +
  theme_bw() + labs(fill = "Predicted\nClass") +
  theme(text = element_text(size=12), axis.text.x = element_text(angle = 30, hjust = 1), legend.position="bottom", legend.direction = "horizontal") +
  scale_fill_manual(values = c("#33A02C", "#B2DF8A" ,"#A6CEE3" ,"#FDBF6F", "#E31A1C" ,"#FB9A99", "#1F78B4")) +
  guides(fill=guide_legend(nrow=2))

grid.arrange(p1,p2, nrow = 1)




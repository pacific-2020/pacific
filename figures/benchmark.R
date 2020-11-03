library(Rmisc)
library(ggplot2)
library(gridExtra)
library(tidyverse)
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

predictions <- read.table("pacific/figures/benchmarkresults.txt", sep = "\t", header = T)
predictionmetrics <- data.frame()
for (did in paste("T",seq(1:100),sep="")){
  for (givenlabel in levels(predictions$knownlabel)) {
    a <- predictions %>% mutate_if(is.factor, as.character) %>% filter( dataid == did )
    predictionmetrics <- bind_rows(predictionmetrics, calcmetrics(a = a, datatype = "all", givenlabel = givenlabel))
    b <- predictions %>% mutate_if(is.factor, as.character) %>% filter( dataid == did & readtype != "EQUAL")
    predictionmetrics <- bind_rows(predictionmetrics, calcmetrics(a = b, datatype = "mm", givenlabel = givenlabel))
    c <- predictions %>% mutate_if(is.factor, as.character) %>% filter( dataid == did & readtype == "EQUAL")
    predictionmetrics <- bind_rows(predictionmetrics, calcmetrics(a = c, datatype = "equal", givenlabel = givenlabel))
  }
}

dp <- read.table("doublepredictions.txt", sep = "\t")
colnames(dp) <- colnames(predictions)
#dp %>% mutate(knownlabel = case_when(knownlabel=="Cornidoviridae" ~ as.factor("Coronaviridae")))
dpm <- data.frame()
for (did in paste("T",seq(1:100),sep="")){
  for (givenlabel in levels(dp$knownlabel)) {
    a <- dp %>% mutate_if(is.factor, as.character) %>% filter( dataid == did )
    dpm <- bind_rows(dpm, calcmetrics(a = a, datatype = "all", givenlabel = givenlabel))
    b <- dp %>% mutate_if(is.factor, as.character) %>% filter( dataid == did & readtype != "EQUAL")
    dpm <- bind_rows(dpm, calcmetrics(a = b, datatype = "mm", givenlabel = givenlabel))
    c <- dp %>% mutate_if(is.factor, as.character) %>% filter( dataid == did & readtype == "EQUAL")
    dpm <- bind_rows(dpm, calcmetrics(a = c, datatype = "equal", givenlabel = givenlabel))
  }
}

###FNR increase
left_join(predictionmetrics %>% filter(datatype=="all") %>% select(dataid,class, fnr) %>% summarySE(groupvars = "class", measurevar = "fnr"), dpm %>% filter(datatype=="all") %>% select(dataid,class, fnr) %>% summarySE(groupvars = "class", measurevar = "fnr"), by = c("class")) %>% mutate(ratio = fnr.y/fnr.x)
##FPR decrease
left_join(predictionmetrics %>% filter(datatype=="all") %>% select(dataid,class, fpr) %>% summarySE(groupvars = "class", measurevar = "fpr"), dpm %>% filter(datatype=="all") %>% select(dataid,class, fpr) %>% summarySE(groupvars = "class", measurevar = "fpr"), by = c("class")) %>% mutate(ratio = fpr.x/fpr.y)


##table for precision, recall, and accuracy, fpr, and fnr

dpmtable <- data.frame()

dpmtable <- bind_rows(dpmtable, data.frame(left_join(filter(dpm, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(precision), max = max(precision)),
                                                     filter(dpm, datatype == "all") %>% summarySE(measurevar = "precision", groupvars = c("class"))
)) %>% select(class,precision,ci) %>% mutate_if(is.numeric, ~format(., nsmall = 6)) %>% mutate(test = "precision", op="(",cp=")",s=" ") %>% unite(val, precision, s, op,ci,cp, sep = ""))

dpmtable <- bind_rows(dpmtable, data.frame(left_join(filter(dpm, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(recall), max = max(recall)),
                                                     filter(dpm, datatype == "all") %>% summarySE(measurevar = "recall", groupvars = c("class"))
)) %>% select(class,recall,ci) %>% mutate_if(is.numeric, ~format(., nsmall = 4)) %>% mutate(test = "recall", op="(",cp=")",s=" ") %>% unite(val, recall, s, op,ci,cp, sep = ""))

dpmtable <- bind_rows(dpmtable, data.frame(left_join(filter(dpm, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(accuracy), max = max(accuracy)),
                                                     filter(dpm, datatype == "all") %>% summarySE(measurevar = "accuracy", groupvars = c("class"))
)) %>% select(class,accuracy,ci) %>% mutate_if(is.numeric, ~format(., nsmall = 4)) %>% mutate(test = "accuracy", op="(",cp=")",s=" ") %>% unite(val, accuracy, s, op,ci,cp, sep = ""))

dpmtable <- bind_rows(dpmtable, data.frame(left_join(filter(dpm, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(fpr), max = max(fpr)),
                                                     filter(dpm, datatype == "all") %>% summarySE(measurevar = "fpr", groupvars = c("class"))
)) %>% select(class,fpr,ci) %>% mutate_if(is.numeric, ~format(., nsmall = 3)) %>% mutate(test = "fpr", op="(",cp=")",s=" ") %>% unite(val, fpr, s, op,ci,cp, sep = ""))


dpmtable <- bind_rows(dpmtable, data.frame(left_join(filter(dpm, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(fnr), max = max(fnr)),
                                                     filter(dpm, datatype == "all") %>% summarySE(measurevar = "fnr", groupvars = c("class"))
)) %>% select(class,fnr,ci) %>% mutate_if(is.numeric, ~format(., nsmall = 1)) %>% mutate(test = "fnr", op="(",cp=")",s=" ") %>% unite(val, fnr, s, op,ci,cp, sep = ""))

spread(dpmtable, test, val)
write.table(spread(dpmtable, test, val), file="performance.metrics.dp.txt", sep = "\t", col.names = T, row.names = F, quote = F)


##table for precision, recall, and accuracy

##table for precision, recall, and accuracy
metricstable <- data.frame()

metricstable <- bind_rows(metricstable, data.frame(left_join(filter(predictionmetrics, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(precision), max = max(precision)),
                     filter(predictionmetrics, datatype == "all") %>% summarySE(measurevar = "precision", groupvars = c("class"))
)) %>% select(class,precision,ci) %>% mutate_if(is.numeric, ~format(., nsmall = 1)) %>% mutate(test = "precision", op="(",cp=")",s=" ") %>% unite(val, precision, s, op,ci,cp, sep = ""))

metricstable <- bind_rows(metricstable, data.frame(left_join(filter(predictionmetrics, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(recall), max = max(recall)),
                     filter(predictionmetrics, datatype == "all") %>% summarySE(measurevar = "recall", groupvars = c("class"))
)) %>% select(class,recall,ci) %>% mutate_if(is.numeric, ~format(., nsmall = 1)) %>% mutate(test = "recall", op="(",cp=")",s=" ") %>% unite(val, recall, s, op,ci,cp, sep = ""))

metricstable <- bind_rows(metricstable, data.frame(left_join(filter(predictionmetrics, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(accuracy), max = max(accuracy)),
                     filter(predictionmetrics, datatype == "all") %>% summarySE(measurevar = "accuracy", groupvars = c("class"))
)) %>% select(class,accuracy,ci) %>% mutate_if(is.numeric, ~format(., nsmall = 4)) %>% mutate(test = "accuracy", op="(",cp=")",s=" ") %>% unite(val, accuracy, s, op,ci,cp, sep = ""))

metricstable <- bind_rows(metricstable, data.frame(left_join(filter(predictionmetrics, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(fpr), max = max(fpr)),
                     filter(predictionmetrics, datatype == "all") %>% summarySE(measurevar = "fpr", groupvars = c("class"))
)) %>% select(class,fpr,ci) %>% mutate_if(is.numeric, ~format(., nsmall = 1)) %>% mutate(test = "fpr", op="(",cp=")",s=" ") %>% unite(val, fpr, s, op,ci,cp, sep = ""))


metricstable <- bind_rows(metricstable, data.frame(left_join(filter(predictionmetrics, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(fnr), max = max(fnr)),
                     filter(predictionmetrics, datatype == "all") %>% summarySE(measurevar = "fnr", groupvars = c("class"))
)) %>% select(class,fnr,ci) %>% mutate_if(is.numeric, ~format(., nsmall = 1)) %>% mutate(test = "fnr", op="(",cp=")",s=" ") %>% unite(val, fnr, s, op,ci,cp, sep = ""))


spread(metricstable, test, val)

write.table(spread(metricstable, test, val), file="performance.metrics.txt", sep = "\t", col.names = T, row.names = F, quote = F)


data.frame(left_join(filter(predictionmetrics, datatype == "mm") %>% group_by(class) %>% dplyr::summarise(min = min(precision), max = max(precision)),
                     filter(predictionmetrics, datatype == "mm") %>% summarySE(measurevar = "precision", groupvars = c("class"))
))

data.frame(left_join(filter(predictionmetrics, datatype == "mm") %>% group_by(class) %>% dplyr::summarise(min = min(recall), max = max(recall)),
                     filter(predictionmetrics, datatype == "mm") %>% summarySE(measurevar = "recall", groupvars = c("class"))
))

data.frame(left_join(filter(predictionmetrics, datatype == "mm") %>% group_by(class) %>% dplyr::summarise(min = min(accuracy), max = max(accuracy)),
                     filter(predictionmetrics, datatype == "mm") %>% summarySE(measurevar = "accuracy", groupvars = c("class"))
))

##fpr and fnr numbers
##heatmap for false positive numbers
predictions %>% 
  filter(as.character(knownlabel) != as.character(predictedlabel)) %>% 
  group_by(dataid,knownlabel,predictedlabel) %>% 
  dplyr::summarise(totalcounts = sum(readcounts)) %>% 
  group_by(dataid, predictedlabel) %>% 
  dplyr::summarise(knownlabel=knownlabel, fppid = totalcounts*100/sum(totalcounts), fpclass = totalcounts,  totalfp= sum(totalcounts)) %>% 
  ungroup() %>% 
  select(predictedlabel, knownlabel, fppid) %>% 
  group_by(predictedlabel, knownlabel) %>% 
  dplyr::summarise(m=mean(fppid, na.rm = T)) %>% 
  filter(! predictedlabel %in% c("pDiscarded", "Discarded")) %>%
  ggplot(aes(x=knownlabel,y=predictedlabel, fill = m)) + 
  geom_tile(aes(fill=m)) + 
  geom_text(aes(label=round(m,1))) +
  scale_fill_gradientn( colours = brewer.pal(n=10, name = "RdYlBu")) +
  theme_bw() +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(fill = "FP %") +
  ylab("Predicted FP Class") +
  xlab("FP Source Class")

###false positive sources bar chart
mycolors <- sample(brewer.pal(n=12,name="Paired"), size = 7, replace = F )
print(mycolors)
p1 <- dp %>% 
  filter(as.character(knownlabel) != as.character(predictedlabel)) %>% 
  group_by(dataid,knownlabel,predictedlabel) %>% 
  dplyr::summarise(totalcounts = sum(readcounts)) %>% 
  group_by(knownlabel, predictedlabel) %>% 
  dplyr::summarise(average = sum(totalcounts)/100) %>% 
  filter(! predictedlabel %in% c("Discarded", "Human", "pDiscarded", "rc_discarded")) %>% 
  ggplot(aes(x=predictedlabel, y = average, fill = knownlabel)) + 
  geom_bar(stat="identity", position = "stack") + 
  ylab("Mean false positive reads") + 
  xlab("Predicted class") +
  theme_bw() + labs(fill = "True\nClass") +
  theme(text = element_text(size=12), axis.text.x = element_text(angle = 30, hjust = 1), legend.position="bottom") +
  scale_fill_manual(values = c("#33A02C", "#B2DF8A" ,"#A6CEE3" ,"#FDBF6F", "#E31A1C" ,"#FB9A99", "#1F78B4"))

p2 <- dp %>% 
  filter(as.character(knownlabel) != as.character(predictedlabel)) %>% 
  group_by(dataid,knownlabel,predictedlabel) %>% 
  dplyr::summarise(totalcounts = sum(readcounts)) %>% 
  group_by(knownlabel, predictedlabel) %>% 
  dplyr::summarise(average = sum(totalcounts)/100) %>% 
  filter(predictedlabel != "Discarded" & average > 0 & knownlabel != "Unrelated" & knownlabel != "Human") %>% 
  ggplot(aes(x=knownlabel, y = average, fill = predictedlabel)) + 
  geom_bar(stat="identity", position = "stack") + 
  ylab("Mean false negative reads") + 
  xlab("True class") +
  theme_bw() + labs(fill = "Predicted\nClass") +
  theme(text = element_text(size=12), axis.text.x = element_text(angle = 30, hjust = 1), legend.position="bottom", legend.direction = "horizontal") +
  scale_fill_manual(values = c("#33A02C", "#B2DF8A" ,"#A6CEE3" ,"#FDBF6F", "#E31A1C" ,"#FB9A99", "#1F78B4")) +
  guides(fill=guide_legend(nrow=2))

grid.arrange(p1,p2, nrow = 1)

###false negative mismatch read percentages

fnmm <- dp %>% filter(as.character(knownlabel) != as.character(predictedlabel) & 
                        predictedlabel != "Discarded" & knownlabel != "Unrelated") %>% 
  mutate(readtype = case_when(readtype != "EQUAL" ~ "MM", TRUE ~ as.character(readtype)))  %>% 
  group_by(knownlabel, readtype) %>% 
  dplyr::summarise(t = sum(readcounts)) %>% ungroup()

left_join(fnmm, fnmm %>% group_by(knownlabel) %>% dplyr::summarise(s = sum(t)), by = c("knownlabel")) %>% mutate(p = t/s) %>%
  ggplot(aes(x=knownlabel,y=p,fill=readtype)) + geom_bar(stat="identity", position = "stack")


dp %>% 
  filter(as.character(knownlabel) != as.character(predictedlabel)) %>% mutate(readtype = case_when(readtype != "EQUAL" ~ "MM", TRUE ~ as.character(readtype))) %>%
  group_by(dataid,knownlabel,predictedlabel, readtype) %>% 
  dplyr::summarise(totalcounts = sum(readcounts)) %>% 
  group_by(knownlabel, predictedlabel, readtype) %>% 
  dplyr::summarise(average = sum(totalcounts)/100) %>% 
  filter(predictedlabel != "Discarded" & average > 0 & knownlabel != "Unrelated") %>%
  ggplot(aes(x=knownlabel, y = average, fill = readtype)) + 
  geom_bar(stat="identity", position = "stack") + 
  ylab("Mean false negative reads") + 
  xlab("Known class") +
  theme_bw() + labs(fill = "Read\ntype") +
  theme(text = element_text(size=12), axis.text.x = element_text(angle = 30, hjust = 1), legend.position="bottom", legend.direction = "horizontal") +
  scale_fill_manual(values = c("#33A02C", "#B2DF8A" ,"#A6CEE3" ,"#FDBF6F", "#E31A1C" ,"#FB9A99", "#1F78B4")) +
  scale_y_continuous(labels = scales::percent)



fprplot <- bind_rows(predictionmetrics %>% 
            filter(datatype=="all" & ! class %in% c("Human", "Unrelated")) %>% 
            select(dataid, class, fpr) %>% mutate(ptype = "p1"), 
          dpm %>% filter(datatype=="all" & ! class %in% c("Human", "Unrelated")) %>% 
            select(dataid, class,fpr) %>% mutate(ptype = "p2")) %>% 
  gather("test", "value", -dataid, -class, -ptype) %>%
  mutate(test = case_when(test == "fpr" ~ "False positive rate", TRUE ~ "nothing")) %>%
  ggplot(aes(x = class, y = value, fill = ptype)) + 
  geom_boxplot(notch = T) + 
  labs(fill = "Prediction\ntype") +
  scale_fill_discrete(labels = c("Normal", "Double")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), text = element_text(size=16)) +
  ylab("False positive rate") +
  xlab("Virus class") +
  ggtitle("False positive rate")

fnrplot <- bind_rows(predictionmetrics %>% 
                       filter(datatype=="all" & ! class %in% c("Human", "Unrelated")) %>% 
                       select(dataid, class, fnr) %>% mutate(ptype = "p1"), 
                     dpm %>% filter(datatype=="all" & ! class %in% c("Human", "Unrelated")) %>% 
                       select(dataid, class,fnr) %>% mutate(ptype = "p2")) %>% 
  gather("test", "value", -dataid, -class, -ptype) %>%
  mutate(test = case_when(test == "fnr" ~ "False negative rate", TRUE ~ "nothing")) %>%
  ggplot(aes(x = class, y = value, fill = ptype)) + 
  geom_boxplot(notch = T) + 
  labs(fill = "Prediction\ntype") +
  scale_fill_discrete(labels = c("Normal", "Double")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), text = element_text(size=16)) +
  ylab("False negative rate") +
  xlab("Virus class") +
  ggtitle("False negative rate")

precisionplot <- bind_rows(predictionmetrics %>% 
                             filter(datatype=="all" & ! class %in% c("Human", "Unrelated")) %>% 
                             select(dataid, class, precision) %>% mutate(ptype = "p1"), 
                           dpm %>% filter(datatype=="all" & ! class %in% c("Human", "Unrelated")) %>% 
                             select(dataid, class,precision) %>% mutate(ptype = "p2")) %>% 
  gather("test", "value", -dataid, -class, -ptype) %>%
  mutate(test = case_when(test == "precision" ~ "Precision", TRUE ~ "nothing")) %>%
  ggplot(aes(x = class, y = value, fill = ptype)) + 
  geom_boxplot(notch = T) + 
  labs(fill = "Prediction\ntype") +
  scale_fill_discrete(labels = c("Normal", "Double")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), text = element_text(size=16)) +
  ylab("Precision") +
  xlab("Virus class") +
  ggtitle("Precision")

accuracyplot <- bind_rows(predictionmetrics %>% 
                             filter(datatype=="all" & ! class %in% c("Human", "Unrelated")) %>% 
                             select(dataid, class, accuracy) %>% mutate(ptype = "p1"), 
                           dpm %>% filter(datatype=="all" & ! class %in% c("Human", "Unrelated")) %>% 
                             select(dataid, class,accuracy) %>% mutate(ptype = "p2")) %>% 
  gather("test", "value", -dataid, -class, -ptype) %>%
  mutate(test = case_when(test == "accuracy" ~ "Accuracy", TRUE ~ "nothing")) %>%
  ggplot(aes(x = class, y = value, fill = ptype)) + 
  geom_boxplot(notch = T) + 
  labs(fill = "Prediction\ntype") +
  scale_fill_discrete(labels = c("Normal", "Double")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), text = element_text(size=16)) +
  ylab("Accuracy") +
  xlab("Virus class") +
  ggtitle("Accuracy") +
  scale_y_continuous(labels = comma)

recallplot <- bind_rows(predictionmetrics %>% 
                          filter(datatype=="all" & ! class %in% c("Human", "Unrelated")) %>% 
                          select(dataid, class, recall) %>% mutate(ptype = "p1"), 
                        dpm %>% filter(datatype=="all" & ! class %in% c("Human", "Unrelated")) %>% 
                          select(dataid, class,recall) %>% mutate(ptype = "p2")) %>% 
  gather("test", "value", -dataid, -class, -ptype) %>%
  mutate(test = case_when(test == "recall" ~ "Recall", TRUE ~ "nothing")) %>%
  ggplot(aes(x = class, y = value, fill = ptype)) + 
  geom_boxplot(notch = T) + 
  labs(fill = "Prediction\ntype") +
  scale_fill_discrete(labels = c("Normal", "Double")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), text = element_text(size=16)) +
  ylab("Recall") +
  xlab("Virus class") +
  ggtitle("Recall")


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



grid_arrange_shared_legend(fprplot,fnrplot,accuracyplot,precisionplot,recallplot, nrow = 2, ncol = 3)


predictionmetrics %>% filter(datatype == "all") %>% mutate(pdp = pdiscarded*100/total, dp = discarded*100/total, tdp = (discarded + pdiscarded)*100/total) %>% summarySE(measurevar = "tdp", groupvars = c("class"))
predictionmetrics %>% filter(datatype == "all") %>% mutate(pdp = pdiscarded*100/total, dp = discarded*100/total, tdp = (discarded + pdiscarded)*100/total) %>% summarySE(measurevar = "dp", groupvars = c("class"))
predictionmetrics %>% filter(datatype == "all") %>% mutate(pdp = pdiscarded*100/total, dp = discarded*100/total, tdp = (discarded + pdiscarded)*100/total) %>% summarySE(measurevar = "pdp", groupvars = c("class"))


predictionmetrics %>% filter(datatype == "mm") %>% summarySE(measurevar = "precision", groupvars = c("class")) %>% arrange(precision)
predictionmetrics %>% filter(datatype == "mm") %>% summarySE(measurevar = "recall", groupvars = c("class")) %>% arrange(recall)
predictionmetrics %>% filter(datatype == "mm") %>% summarySE(measurevar = "accuracy", groupvars = c("class")) %>% arrange(accuracy)





data.frame(left_join(filter(predictionmetrics, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(fpr), max = max(fpr)),
          filter(predictionmetrics, datatype == "all") %>% summarySE(measurevar = "fpr", groupvars = c("class"))
))




          
predictionmetrics %>% filter(datatype == "all") %>% summarySE(measurevar = "fpr", groupvars = c("class")) %>% arrange(fpr)

scaleFUN <- function(x) sprintf("%.4f", x)
metricsplot <- dpm %>% filter(class != "Human" & class != "Unrelated") %>% 
  select(dataid, class, datatype, precision, recall, accuracy) %>%
  gather( key = "test", value = "metric", -dataid, -class, -datatype) %>% 
  mutate(test = case_when(test=="precision" ~ "Precision", test == "accuracy" ~ "Accuracy", test == "recall" ~ "Recall", TRUE ~ as.character(test))) %>%
  mutate(test = as.factor(test)) %>%
  arrange(dataid,class) %>% 
  ggplot(aes(x=class,y=metric,fill=datatype)) + 
  geom_boxplot(notch = T, lwd=0.2) + 
  facet_wrap(.~test, scales = "free_y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  ylab("Value") +
  labs(x = NULL, fill = "Read types") +
  scale_fill_discrete(labels = c("All", "Mismatch", "Exact")) +
  theme(text = element_text(size=20)) + scale_y_continuous(labels=scaleFUN)


fprfnrplot <- dpm %>% filter(class != "Human" & class != "Unrelated") %>% 
  select(dataid, class, datatype, fpr, fnr) %>% 
  gather("metric","value", -dataid, -class, -datatype) %>%
  mutate(metric = case_when(metric=="fpr" ~ "False positive rate", metric == "fnr" ~ "False negative rate", TRUE ~ as.character(metric))) %>%
  mutate(metric = as.factor(metric)) %>%
  ggplot(aes(x=class,y=value,fill=datatype)) + 
  geom_boxplot(notch = T, lwd=0.2) + 
  facet_wrap(.~metric, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  xlab("Virus class") + ylab("Value") + 
  scale_fill_discrete(labels = c("All", "Mismatch", "Exact")) + 
  labs(fill = "Read types") + 
  theme(text = element_text(size=20))

grid.arrange(metricsplot,fprfnrplot, nrow = 2)






predictionmetrics %>% filter(class != "Human" & class != "Unrelated") %>% 
  select(dataid, class, datatype, fpr, fnr) %>% 
  gather
  ggplot(aes(x=class,y=fpr,fill=datatype)) + 
  geom_boxplot(notch = T) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  


##to see which class is introducing FP in which class
predictions %>% 
  filter(predictedlabel != "Discarded" & predictedlabel != "pDiscarded" & as.character(knownlabel) != as.character(predictedlabel) & predictedlabel != "Human") %>% 
  group_by(knownlabel, predictedlabel) %>% 
  dplyr::summarise(s = sum(readcounts)/100) %>% 
  ggplot(aes(x=knownlabel, y = predictedlabel, fill = s)) + 
  geom_tile() + 
  scale_fill_viridis_b() + 
  theme_bw()


predictions %>% filter(as.character(knownlabel) != as.character(predictedlabel)) %>% 
  group_by(dataid,knownlabel,predictedlabel) %>% 
  dplyr::summarise(totalcounts = sum(readcounts)) %>% 
  group_by(dataid, predictedlabel) %>% 
  dplyr::summarise( knownlabel=knownlabel, fppid = totalcounts*100/sum(totalcounts), fpclass = totalcounts,  totalfp= sum(totalcounts)) %>% ggplot(aes(x=knownlabel,y=predictedlabel, fill=fppid)) + geom_tile() + scale_fill_viridis_c()


predictions %>% 
  filter(as.character(knownlabel) != as.character(predictedlabel)) %>% 
  group_by(dataid,knownlabel,predictedlabel) %>% 
  dplyr::summarise(totalcounts = sum(readcounts)) %>% 
  group_by(dataid, predictedlabel) %>% 
  dplyr::summarise(knownlabel=knownlabel, fppid = totalcounts*100/sum(totalcounts), fpclass = totalcounts,  totalfp= sum(totalcounts)) %>% 
  ggplot(aes(x=knownlabel,y=predictedlabel, fill=fppid)) + 
  geom_tile() + 
  stat_summary(fun.data = fun_mean, geom="tile") + 
  scale_fill_gradientn( colours = brewer.pal(n=10, name = "RdYlBu")) +
  theme_bw()

#  mutate(fppid = replace_na(fppid, 0)) %>%



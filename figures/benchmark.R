library(Rmisc)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)


calcmetrics <- function(a, datatype, givenlabel) {

  ##metrics calculated according to definitions at https://en.wikipedia.org/wiki/Precision_and_recall
  pdiscarded <- sum(a %>% filter(knownlabel == givenlabel & predictedlabel == "pDiscarded") %>% .$readcounts)
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
                       discarded = discarded,
                       pdiscarded = pdiscarded
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


##table for precision, recall, and accuracy

##table for precision, recall, and accuracy

data.frame(left_join(filter(predictionmetrics, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(precision), max = max(precision)),
                     filter(predictionmetrics, datatype == "all") %>% summarySE(measurevar = "precision", groupvars = c("class"))
))

data.frame(left_join(filter(predictionmetrics, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(recall), max = max(recall)),
                     filter(predictionmetrics, datatype == "all") %>% summarySE(measurevar = "recall", groupvars = c("class"))
))

data.frame(left_join(filter(predictionmetrics, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(accuracy), max = max(accuracy)),
                     filter(predictionmetrics, datatype == "all") %>% summarySE(measurevar = "accuracy", groupvars = c("class"))
))

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

data.frame(left_join(filter(predictionmetrics, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(fpr), max = max(fpr)),
                     filter(predictionmetrics, datatype == "all") %>% summarySE(measurevar = "fpr", groupvars = c("class"))
))

data.frame(left_join(filter(predictionmetrics, datatype == "all") %>% group_by(class) %>% dplyr::summarise(min = min(fnr), max = max(fnr)),
                     filter(predictionmetrics, datatype == "all") %>% summarySE(measurevar = "fnr", groupvars = c("class"))
))





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
metricsplot <- predictionmetrics %>% filter(class != "Human" & class != "Unrelated") %>% 
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


fprfnrplot <- predictionmetrics %>% filter(class != "Human" & class != "Unrelated") %>% 
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
  theme_bw()


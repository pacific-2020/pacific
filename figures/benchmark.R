library(Rmisc)
library(ggplot2)
library(gridExtra)
library(tidyverse)

calcmetrics <- function(a, datatype, givenlabel) {

  ##metrics calculated according to definitions at https://en.wikipedia.org/wiki/Precision_and_recall
  pdiscarded <- sum(a %>% filter(knownlabel == givenlabel & predictedlabel == "pDiscarded") %>% .$readcounts)
  discarded  <- sum(a %>% filter(knownlabel == givenlabel & predictedlabel == "Discarded") %>% .$readcounts)
  totalreads <- sum(a %>% filter(knownlabel == givenlabel) %>% .$readcounts)
  a <- a %>% filter(predictedlabel != "Discarded" & predictedlabel != "pDiscarded")
  
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

predictionmetrics %>% filter(datatype == "all") %>% mutate(pdp = pdiscarded*100/total, dp = discarded*100/total, tdp = (discarded + pdiscarded)*100/total) %>% summarySE(measurevar = "tdp", groupvars = c("class"))
predictionmetrics %>% filter(datatype == "all") %>% mutate(pdp = pdiscarded*100/total, dp = discarded*100/total, tdp = (discarded + pdiscarded)*100/total) %>% summarySE(measurevar = "dp", groupvars = c("class"))
predictionmetrics %>% filter(datatype == "all") %>% mutate(pdp = pdiscarded*100/total, dp = discarded*100/total, tdp = (discarded + pdiscarded)*100/total) %>% summarySE(measurevar = "pdp", groupvars = c("class"))

predictionmetrics %>% filter(datatype == "all") %>% summarySE(measurevar = "precision", groupvars = c("class")) %>% arrange(precision)
predictionmetrics %>% filter(datatype == "all") %>% summarySE(measurevar = "recall", groupvars = c("class")) %>% arrange(recall)
predictionmetrics %>% filter(datatype == "all") %>% summarySE(measurevar = "accuracy", groupvars = c("class")) %>% arrange(accuracy)


predictionmetrics %>% filter(datatype == "mm") %>% summarySE(measurevar = "precision", groupvars = c("class")) %>% arrange(precision)
predictionmetrics %>% filter(datatype == "mm") %>% summarySE(measurevar = "recall", groupvars = c("class")) %>% arrange(recall)
predictionmetrics %>% filter(datatype == "mm") %>% summarySE(measurevar = "accuracy", groupvars = c("class")) %>% arrange(accuracy)

left_join(filter(predictionmetrics, datatype == "mm") %>% group_by(class) %>% dplyr::summarise(min = min(precision), max = max(precision)),
          filter(predictionmetrics, datatype == "mm") %>% summarySE(measurevar = "precision", groupvars = c("class"))
          )

left_join(filter(predictionmetrics, datatype == "mm") %>% group_by(class) %>% dplyr::summarise(min = min(balancedaccuracy), max = max(balancedaccuracy)),
          filter(predictionmetrics, datatype == "mm") %>% summarySE(measurevar = "balancedaccuracy", groupvars = c("class"))
)





          
predictionmetrics %>% filter(datatype == "all") %>% summarySE(measurevar = "fpr", groupvars = c("class")) %>% arrange(fpr)


predictionmetrics %>% filter(class != "Human" & class != "Unrelated") %>% 
  select(dataid, class, datatype, precision, recall, accuracy) %>% 
  gather( key = "test", value = "metric", -dataid, -class, -datatype) %>% 
  arrange(dataid,class) %>% 
  ggplot(aes(x=class,y=metric,fill=datatype)) + 
  geom_boxplot(notch = T) + 
  facet_wrap(.~test) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

predictionmetrics %>% filter(class != "Human" & class != "Unrelated") %>% 
  select(dataid, class, datatype, fpr) %>% 
  ggplot(aes(x=class,y=fpr,fill=datatype)) + 
  geom_boxplot(notch = T) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


##to see which class is introducing FP in which class
predictions %>% 
  filter(predictedlabel != "Discarded" & predictedlabel != "pDiscarded" & as.character(knownlabel) != as.character(predictedlabel) & predictedlabel != "Human") %>% 
  group_by(knownlabel, predictedlabel) %>% 
  dplyr::summarise(s = sum(readcounts)/100) %>% 
  ggplot(aes(x=knownlabel, y = predictedlabel, fill = s)) + 
  geom_tile() + 
  scale_fill_viridis_b() + 
  theme_bw()









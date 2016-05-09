#####################################
### BISC577A Unit3 - Assignment 2 ###
###           Shuo Li             ###
#####################################

setwd("/Users/lishuo/Desktop/577")

### (3) preparation of high-throughput in vitra data analysis

## install Bioconductor
#source("https://bioconductor.org/biocLite.R")
#biocLite()

## install DNAshapeR
#biocLite("DNAshapeR")

## install caret
#install.packages("caret")
library(DNAshapeR)
library(caret)

## import data
mad.fasta <- "Mad.txt.fa.txt"
mad.pred <- getShape(mad.fasta)
max.fasta <- "Max.txt.fa.txt"
max.pred <- getShape(max.fasta)
myc.fasta <- "Myc.txt.fa.txt"
myc.pred <- getShape(myc.fasta)




### (4) build prediction models for in vitro data
f1mer <- c("1-mer")
f1mershape <- c("1-mer","1-shape")

## Mad
# generate feature vector
mad.1mer <- encodeSeqShape(mad.fasta, mad.pred, f1mer)
mad.1mershape <- encodeSeqShape(mad.fasta, mad.pred, f1mershape)
# build MLR model with Caret
mad.exp <- "Mad.txt"
mad.expdata <- read.table(mad.exp)
mad.df1mer <- data.frame(affinity = mad.expdata$V2, mad.1mer)
mad.df1mershape <- data.frame(affinity = mad.expdata$V2, mad.1mershape)
# arguments setting for Caret
mad.train <- trainControl(method = "cv", number = 10, savePredictions = T)
# prediction with L2-regularization
mad.model1mer <- train(affinity~., data = mad.df1mer, trControl=mad.train, 
                 method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
mad.model1mershape <- train(affinity~., data = mad.df1mershape, trControl=mad.train, 
              method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
# plot R^2 curve
jpeg("Mad.jpeg")
plot(x = c(-15:15), y = mad.model1mer$results$Rsquared, ylim = c(0.4,0.9), xlab = "log2(lambda)", ylab = "R^2",
     main = "Mad prediction with L2-regularization", type = 'l', col = 4)
lines(x = c(-15:15), y = mad.model1mershape$results$Rsquared, col = 2)
legend("topright", legend = c("1mer        ","1mer+shape        "), lty = c(1,1), col = c(4,2))
dev.off()

## Max
# generate feature vector
max.1mer <- encodeSeqShape(max.fasta, max.pred, f1mer)
max.1mershape <- encodeSeqShape(max.fasta, max.pred, f1mershape)
# build MLR model with Caret
max.exp <- "Max.txt"
max.expdata <- read.table(max.exp)
max.df1mer <- data.frame(affinity = max.expdata$V2, max.1mer)
max.df1mershape <- data.frame(affinity = max.expdata$V2, max.1mershape)
# arguments setting for Caret
max.train <- trainControl(method = "cv", number = 10, savePredictions = T)
# prediction with L2-regularization
max.model1mer <- train(affinity~., data = max.df1mer, trControl=max.train, 
                       method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
max.model1mershape <- train(affinity~., data = max.df1mershape, trControl=max.train, 
                            method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
# plot R^2 curve
jpeg("Max.jpeg")
plot(x = c(-15:15), y = max.model1mer$results$Rsquared, ylim = c(0.4,0.9), xlab = "log2(lambda)", ylab = "R^2",
     main = "Max prediction with L2-regularization", type = 'l', col = 4)
lines(x = c(-15:15), y = max.model1mershape$results$Rsquared, col = 2)
legend("topright", legend = c("1mer         ","1mer+shape         "), lty = c(1,1), col = c(4,2))
dev.off()

## Myc
# generate feature vector
myc.1mer <- encodeSeqShape(myc.fasta, myc.pred, f1mer)
myc.1mershape <- encodeSeqShape(myc.fasta, myc.pred, f1mershape)
# build MLR model with Caret
myc.exp <- "Myc.txt"
myc.expdata <- read.table(myc.exp)
myc.df1mer <- data.frame(affinity = myc.expdata$V2, myc.1mer)
myc.df1mershape <- data.frame(affinity = myc.expdata$V2, myc.1mershape)
# arguments setting for Caret
myc.train <- trainControl(method = "cv", number = 10, savePredictions = T)
# prediction with L2-regularization
myc.model1mer <- train(affinity~., data = myc.df1mer, trControl=myc.train, 
                       method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
myc.model1mershape <- train(affinity~., data = myc.df1mershape, trControl=myc.train, 
                            method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
# plot R^2 curve
jpeg("Myc.jpeg")
plot(x = c(-15:15), y = myc.model1mer$results$Rsquared, ylim = c(0.4,0.9), xlab = "log2(lambda)", ylab = "R^2",
     main = "Myc prediction with L2-regularization", type = 'l', col = 4)
lines(x = c(-15:15), y = myc.model1mershape$results$Rsquared, col = 2)
legend("topright", legend = c("1mer        ","1mer+shape        "), lty = c(1,1), col = c(4,2))
dev.off()



### (5) high-throughput in vitro data analysis
R2.1mer <- c(mad.model1mer$result$Rsquared[which(mad.model1mer$results$lambda == mad.model1mer$bestTune$lambda)],
             max.model1mer$result$Rsquared[which(max.model1mer$results$lambda == max.model1mer$bestTune$lambda)],
             myc.model1mer$result$Rsquared[which(myc.model1mer$results$lambda == myc.model1mer$bestTune$lambda)])
R2.1mershape <- c(mad.model1mershape$result$Rsquared[which(mad.model1mershape$results$lambda == mad.model1mershape$bestTune$lambda)],
             max.model1mershape$result$Rsquared[which(max.model1mershape$results$lambda == max.model1mershape$bestTune$lambda)],
             myc.model1mershape$result$Rsquared[which(myc.model1mershape$results$lambda == myc.model1mershape$bestTune$lambda)])

## performance comparison
jpeg("performance.jpeg")
plot(x = R2.1mer, y = R2.1mershape, xlim = c(0.75,0.9), 
     ylim = c(0.75,0.9), pch = c(17,18,20), col = c(2,3,4),
     xlab = "R^2 1mer", ylab = "R^2 1mer+shape",
     main = "preformance of 1mer and 1mer+shape")
lines(x = c(0.75,0.9), y = c(0.75,0.9), lty = 2, col = "dark gray")
legend("bottomright", legend = c("Mad    ", "Max    ", "Myc    "),
       col = c(2,3,4), pch = c(17,18,20))
dev.off()

## compute p-value
wilcox.test(x = R2.1mer, y = R2.1mershape, alternative = "less") #alternative hypothesis: x is shifted to left of y


### (6) preparation of high-throughput in vivo data analysis

## install AnnotationHub
#source("https://bioconductor.org/biocLite.R")
#biocLite("AnnotationHub")
#biocLite("BSgenome.Mmusculus.UCSC.mm10")
library("AnnotationHub")

## find CTCF data
library(BSgenome.Mmusculus.UCSC.mm10)
ah <- AnnotationHub()
CTCF <- ah[["AH28451"]]

## generate fasta file
seqlevelsStyle(CTCF) <- "UCSC"
getFasta(GR = CTCF, BSgenome = Mmusculus, width = 401, filename = "CTCF.fa")


### (7) High-throughput in vivo data analysis

## ensemble plots
CTCFshape <- getShape("CTCF.fa")
jpeg("plotMGW.jpeg")
plotShape(CTCFshape$MGW,main = "plotShape of MGW")
dev.off()
jpeg("plotProT.jpeg")
plotShape(CTCFshape$ProT,main = "plotShape of ProT")
dev.off()
jpeg("plotRoll.jpeg")
plotShape(CTCFshape$Roll,main = "plotShape of Roll")
dev.off()
jpeg("plotHelT.jpeg")
plotShape(CTCFshape$HelT,main = "plotShape of HelT")
dev.off()
jpeg("heatMGW.jpeg")
heatShape(CTCFshape$MGW, 401,main = "heatShape of MGW")
dev.off()
jpeg("heatProT.jpeg")
heatShape(CTCFshape$ProT,401,main = "heatShape of ProT")
dev.off()
jpeg("heatRoll.jpeg")
heatShape(CTCFshape$Roll,400,main = "heatShape of Roll")
dev.off()
jpeg("heatHelt.jpeg")
heatShape(CTCFshape$HelT,400,main = "heatShape of HelT")
dev.off()



### (8) High-throughput in vivo data analysis

## generate random sequences
chrName <- names(Mmusculus)[1:22]
chrLength <- seqlengths(Mmusculus)[1:22]
randomseq <- GRanges()
while(length(randomseq) < 1000){
  tmpChrName <- sample(chrName, 1)
  tmpChrLength <- chrLength[tmpChrName]
  tmpseq <- GRanges(seqnames = tmpChrName, strand = "+", IRanges(start = sample(1:(tmpChrLength-30),1), width = 30))
  if(length(findOverlaps(CTCF, tmpseq)) == 0) randomseq <- c(randomseq, tmpseq)
}

## trim CTCF sequence
CTCFtrimmed <- GRanges()
start <- sample(1:(length(CTCF)-1000),1)
for(i in start:(start+1000)){ #I trimmed the 1000 consecutive sequences with a random start position
  tmpseq <- GRanges(seqnames = CTCF[i]@seqnames@values, strand = CTCF[i]@strand@values,
                    IRanges(start = CTCF[i]$feature_start_position, width = 30))
  CTCFtrimmed <- c(CTCFtrimmed, tmpseq)
}
#CTCFts <- CTCFtrimmed[sample(1:length(CTCFtrimmed),1000)]

## overlap checking
findOverlaps(CTCFtrimmed, randomseq)

## generate fasta file for trimed CTCF and random sequence
getFasta(GR = CTCFtrimmed, BSgenome = Mmusculus, width = 30, filename = "CTCFtrimmed.fa")
getFasta(GR = randomseq, BSgenome = Mmusculus, width = 30, filename = "randomseq.fa")

## combine two datasets
boundFasta <- readBStringSet("CTCFtrimmed.fa")
unboundFasta <- readBStringSet("randomseq.fa")
names(unboundFasta) <- paste0( names(unboundFasta), "_unbound")
writeXStringSet( c(boundFasta, unboundFasta), "ctcf.fa" )

## generate binding classfication file
boundTxt <- cbind( sapply(1:length(boundFasta), 
                          function(x) as.character(boundFasta[[x]])), 
                   matrix(1, length(boundFasta), 1))
unboundTxt <- cbind( sapply(1:length(unboundFasta),
                            function(x) as.character(unboundFasta[[x]])),
                     matrix(0, length(unboundFasta), 1))
write.table(rbind(boundTxt, unboundTxt), "ctcf.txt", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

## prepare shape and outcome dataframe
shape <- getShape("ctcf.fa")
ctcf1mer <- encodeSeqShape("ctcf.fa", shape, c("1-mer"))
ctcf1mershape <- encodeSeqShape("ctcf.fa", shape, c("1-mer","1-shape"))
exp_data <- read.table("ctcf.txt")
exp_data$V2 <- ifelse(exp_data$V2 == 1 , "Y", "N")
ctcf.df1mer <- data.frame(outcome = as.factor(exp_data$V2), ctcf1mer)
ctcf.df1mershape <- data.frame(outcome = as.factor(exp_data$V2), ctcf1mershape)

## prediction
trainControl <- trainControl(method = "cv", number = 2, 
                             savePredictions = TRUE, classProbs = TRUE)
model.1mer <- train(outcome~ ., data = ctcf.df1mer, trControl = trainControl,
               method = "glm", family = binomial, metric ="ROC")
summary(model.1mer)
model.1mershape <- train(outcome~ ., data = ctcf.df1mershape, trControl = trainControl,
                    method = "glm", family = binomial, metric ="ROC")
summary(model.1mershape)

## plot ROC curve
#install.packages("ROCR")
library("ROCR")
pred1mer <- prediction(model.1mer$pred$Y, model.1mer$pred$obs )
perf1mer <- performance( pred1mer, "tpr", "fpr" )
pred1mershape <- prediction(model.1mershape$pred$Y, model.1mershape$pred$obs )
perf1mershape <- performance( pred1mershape, "tpr", "fpr" )
#summary(performance)
jpeg("ROC.jpeg")
plot(perf1mer,main = "ROC curve of 1-mer model", col = "blue")
lines(x = perf1mershape@x.values[[1]], y = perf1mershape@y.values[[1]], main = "ROC curve of 1-mer+1-shape model", col = "red")
lines(c(0,1),c(0,1),col = "dark gray" , lty =2)
legend("bottomright", legend = c("1mer        ", "1mer+1shape         "),
       col = c("blue","red"), lty = c(1,1))
dev.off()

## caluculate AUC(area under ROC curve)
auc1mer <- performance(pred1mer, "auc")
auc1mer <- unlist(slot(auc1mer, "y.values"))
auc1mershape <- performance(pred1mershape, "auc")
auc1mershape <- unlist(slot(auc1mershape, "y.values"))


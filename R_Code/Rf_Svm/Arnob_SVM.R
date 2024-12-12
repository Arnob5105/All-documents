##################################################
rm(list = ls())
setwd("D:\\Bayazid\\ML\\ML_Arnob")
#require(randomForest)
#require(ggplot2)
require(pROC)
require(caret)
require(combinat)
require(e1071)
library(ROCR)

gse1 <- read.csv(file="DataC_N.csv", header = TRUE, sep=",")
#gse1 <- read.csv(file="Data_RF.csv", header = TRUE, sep=",")

dim(gse1)
gse2 <- t(gse1)
gse3 <- data.frame(gse2[-1,])
class <- c(rep("control",59),rep("case",385));class
gene <- c("class",gse1[,1]);gene

gse4 <- cbind(class,gse3)
colnames(gse4) <- gene
rownames(gse4)<- NULL




d <- dim(gse4);d
var= c()
for (i in 1:d[2]-1)
{
  var[i] <- paste0('var',i)
}
var

varname <- c("class",var)

colnames(gse4) <- varname

data <- gse4

data$class <- as.factor(class)
levels(data$class)= c(0,1)
str(data)

data <- lapply(data,as.numeric)
head(str(data))
data$class <- as.factor(class)
str(data)

#data


Hub <- colnames(data)[-1]; Hub

acc = rep(NA, d[2]-1) 
auc <-  rep(NA, d[2]-1)
data[1]
test_idx = sort(sample(1:d[1], round(d[1]*.4,0)));test_idx   # d[1]=126
test = data[test_idx, ]
train= data[-test_idx, ]
str(test[,c(1,3)])
library(ROCR)
#for (i in 1:d[2]-1){
i=15
  training = train[,c(1,i+1)]
  testing = test[,c(1,i+1)]
  svm.model <- svm(class~.,data = training, kernel = 'linear',type = 'C-classification')
  svm.model
  svm.model.pred<-predict(svm.model,new_data=testing ,decision.values=TRUE)
  svm.model.probstest<-attr(svm.model.pred,"decision.values")
  svm.model.classtest<-predict(svm.model,new_data=testing, type="class")
  pred<-prediction(svm.model.probstest,testing$class)
  svm.model.performance_SVM <- performance(pred,"tpr","fpr")
  svm.auc<-performance(pred,"auc")@y.values[[1]];svm.auc
  auc[i] <- svm.auc
  svm.fpr <- performance(pred,"fpr")@y.values[[1]];svm.fpr
  svm.acc<-performance(pred,"acc")@y.values[[1]];svm.acc
  index_fpr <- which(svm.fpr<=0.10); index_fpr
  ind <- index_fpr[length(index_fpr)]
  acc.final <- svm.acc[ind];acc.final
  acc[i] <- acc.final
 #}

AUC <- round(auc,2)
Accuracy <- round(accuracies,2)

result <- data.frame(cbind(gene[-1],AUC,Accuracy)); result 
library(dplyr)
Output <- result %>% filter(  auc >0.90 & Accuracy>.90 )
Output

write.csv(result,"result_datac.csv")
write.csv(Output,"Output_datac.csv")

##################################################
rm(list = ls())
setwd("D:/Arnob/RF")
require(randomForest)
require(ggplot2)
require(pROC)
require(caret)
require(combinat)

gse1 <- read.csv(file="DataC_N.csv", header = TRUE, sep=",")
dim(gse1)
gse2 <- t(gse1)
gse3 <- data.frame(gse2[-1,])
class <- c(rep("control",59),rep("case",385));class
gene <- c("class",gse1[,1]);gene

gse4 <- cbind(class,gse3)
colnames(gse4) <- gene
rownames(gse4)<- NULL
write.csv(gse4,"Data_RF.csv")
d <- dim(gse4);d
var= c()
for (i in 1:d[2]-1)
{
  var[i] <- paste0('var',i)
}
var

varname <- c("class",var)

colnames(gse4) <- varname
#str(gse4)

data <- gse4

Hub <- colnames(data)[-1]; Hub

accuracies = rep(NA, d[2]-1) 
auc <-  rep(NA, d[2]-1)

test_idx = sort(sample(1:d[1], round(d[1]*.4,0)));test_idx   # d[1]=126
test = data[test_idx, ]
train= data[-test_idx, ]
test[,c(1,3)]

for (i in 1:d[2]-1){
  train_data = train[,c(1,i+1)]
  test_data = test[,c(1,i+1)]
  rf_model = randomForest(as.factor(train_data$class)~., data=train_data)
  prediction = predict(rf_model, newdata=test_data, type="response")
  acc = confusionMatrix(as.factor(test_data$class), prediction)$overall[1]
  accuracies[i] = acc
  require(pROC)
  rf.roc<-roc(as.factor(train_data$class),rf_model$votes[,2])
  auc[i] <- auc(rf.roc)
}

AUC <- round(auc,2)
Accuracy <- round(accuracies,2)

result <- data.frame(cbind(gene[-1],AUC,Accuracy)); result 
library(dplyr)
Output <- result %>% filter(  auc >0.85 & Accuracy>.85 )
Output

write.csv(result,"result_randomf_0.85.csv")
write.csv(Output,"Output_randomf_0.85.csv")

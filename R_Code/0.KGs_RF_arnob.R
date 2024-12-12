##==========================================================##
##************* Random Forest ******************************##
#===========================================================##
rm(list = ls())
setwd("E:\\Arnob\\RF") 

## Install required packages
#install.packages("randomForest")
#install.packages("ROCR")
#install.packages("ROSE")
#install.packages("ggplot2")
#install.packages("pROC")
#install.packages("e1071")
#install.packages("caret")
#install.packages("rpart")
#install.packages("gtools")

## Load the required library
require(randomForest)
require(ROCR)
require(ROSE)
require(ggplot2)
require(pROC)
require("e1071")
require("rpart")
require("gtools")
require(ROC)

### Input Train & Test Data
KGsData <- read.csv(file="main.csv", header =TRUE, sep=",")
Hub <- colnames(KGsData)[-1]
Hub
d <- dim(KGsData);d
#set.seed(100)
train_idx = sort(sample(1:d[1], round(d[1]*.6,0)))
data_train = as.data.frame(KGsData[train_idx, ])

data_test1 = as.data.frame(KGsData[-train_idx, ])
data_test2 =  read.csv(file="case_con_50161.csv", header = TRUE, sep=",")

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Test data 1
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#set.seed(100)
rfFit = randomForest(as.factor(data_train$class)~., data=data_train)
rfFit

Pred.RFF_test_1 <- predict(rfFit, data_test1, "prob")
table(predict(rfFit, data_test1), data_test1$class)

PredictionsTest1 = prediction(Pred.RFF_test_1[,2], data_test1$class)
performance(PredictionsTest1,"tpr","fpr")

####Confusion matrix####
mt1=length(PredictionsTest1@cutoffs[[1]]);mt1
pDE1<-rep(0,mt1); TPR1=rep(0,mt1); TNR1=rep(0,mt1); FPR1=rep(0,mt1); FNR1=rep(0,mt1); 
ER1=rep(0,mt1); FDR1=rep(0,mt1);ACCR1=rep(0,mt1);MCCR1=rep(0,mt1);MCRR1=rep(0,mt1);  
for (ii in 1:mt1) 
{
  TPR1[ii]<-PredictionsTest1@tp[[1]][ii]/(PredictionsTest1@tp[[1]][ii]+PredictionsTest1@fn[[1]][ii])
  FPR1[ii]<-PredictionsTest1@fp[[1]][ii]/(PredictionsTest1@tn[[1]][ii]+PredictionsTest1@fp[[1]][ii])
  TNR1[ii]<-PredictionsTest1@tn[[1]][ii]/(PredictionsTest1@tn[[1]][ii]+PredictionsTest1@fp[[1]][ii])
  FNR1[ii]<-PredictionsTest1@fn[[1]][ii]/(PredictionsTest1@fn[[1]][ii]+PredictionsTest1@tp[[1]][ii])
  ACCR1[ii]<-(PredictionsTest1@tp[[1]][ii]+PredictionsTest1@tn[[1]][ii])/(PredictionsTest1@tn[[1]][ii]+PredictionsTest1@fp[[1]][ii]+PredictionsTest1@tp[[1]][ii]+PredictionsTest1@fn[[1]][ii])
  MCCR1[ii]<-((PredictionsTest1@tp[[1]][ii]*PredictionsTest1@tn[[1]][ii])-(PredictionsTest1@fp[[1]][ii]*PredictionsTest1@fn[[1]][ii]))/sqrt((PredictionsTest1@tp[[1]][ii]+PredictionsTest1@fp[[1]][ii])*(PredictionsTest1@tp[[1]][ii]+PredictionsTest1@fn[[1]][ii])*(PredictionsTest1@tn[[1]][ii]+PredictionsTest1@fp[[1]][ii])*(PredictionsTest1@tn[[1]][ii]+PredictionsTest1@fn[[1]][ii]))
  MCRR1[ii]<-(PredictionsTest1@fp[[1]][ii]+PredictionsTest1@fn[[1]][ii])/(PredictionsTest1@tn[[1]][ii]+PredictionsTest1@fp[[1]][ii]+PredictionsTest1@fn[[1]][ii])
  FDR1[ii]<-PredictionsTest1@fp[[1]][ii]/(PredictionsTest1@fp[[1]][ii]+PredictionsTest1@tp[[1]][ii])
}
R1 <- rocdemo.sca(rbinom(40,1,.3), rnorm(40), dxrule.sca, 
                  caseLabel="new case", markerLabel="demo Marker" )
R1@sens<-TPR1; #Sensitivity
R1@spec<-(1-FPR1);#Specificity 
ACC1<-(TPR1+TNR1)/(TPR1+TNR1+FPR1+FNR1)
MCC1<-((TPR1*TNR1)-(FPR1*FNR1))/sqrt((TPR1+FPR1)*(TPR1+FNR1)*(TNR1+FPR1)*(TNR1+FNR1))
MCR1<-(FPR1+FNR1)/(TPR1+TNR1+FPR1+FNR1)
FDR1<-(FPR1)/(TPR1+FPR1)
FNR1=(FNR1)/(FNR1+TPR1)
FDR1 = (FPR1/(TPR1+FPR1))
AUC1<-print(AUC(R1));   # AUC
pAUC1<-print(pAUC(R1, 0.1))   # pAUC upto FPR<=0.1
FPR1
TPR1[[27]]
TNR1[[27]]
FNR1[[27]]
FPR1[[27]]
FDR1[[27]]
ACC1[[27]]
MCC1[[27]]
MCR1[[27]]


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Train Test_2 RF
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Pred.RFF_test_2 <- predict(rfFit, data_test2, "prob")
table(predict(rfFit, data_test2), data_test2$class)

predictionstest = prediction(Pred.RFF_test_2[,2], data_test2$class)
performance(predictionstest,"tpr","fpr")

####Confusion matrix####
m=length(predictionstest@cutoffs[[1]]);m
pDEo<-rep(0,m); TPRo=rep(0,m); TNRo=rep(0,m); FPRo=rep(0,m); FNRo=rep(0,m); 
ERo=rep(0,m); FDRo=rep(0,m);ACCRo=rep(0,m);MCCRo=rep(0,m);MCRRo=rep(0,m);  
for (ii in 1:m) 
{
  TPRo[ii]<-predictionstest@tp[[1]][ii]/(predictionstest@tp[[1]][ii]+predictionstest@fn[[1]][ii])
  FPRo[ii]<-predictionstest@fp[[1]][ii]/(predictionstest@tn[[1]][ii]+predictionstest@fp[[1]][ii])
  TNRo[ii]<-predictionstest@tn[[1]][ii]/(predictionstest@tn[[1]][ii]+predictionstest@fp[[1]][ii])
  FNRo[ii]<-predictionstest@fn[[1]][ii]/(predictionstest@fn[[1]][ii]+predictionstest@tp[[1]][ii])
  ACCRo[ii]<-(predictionstest@tp[[1]][ii]+predictionstest@tn[[1]][ii])/(predictionstest@tn[[1]][ii]+predictionstest@fp[[1]][ii]+predictionstest@tp[[1]][ii]+predictionstest@fn[[1]][ii])
  MCCRo[ii]<-((predictionstest@tp[[1]][ii]*predictionstest@tn[[1]][ii])-(predictionstest@fp[[1]][ii]*predictionstest@fn[[1]][ii]))/sqrt((predictionstest@tp[[1]][ii]+predictionstest@fp[[1]][ii])*(predictionstest@tp[[1]][ii]+predictionstest@fn[[1]][ii])*(predictionstest@tn[[1]][ii]+predictionstest@fp[[1]][ii])*(predictionstest@tn[[1]][ii]+predictionstest@fn[[1]][ii]))
  MCRRo[ii]<-(predictionstest@fp[[1]][ii]+predictionstest@fn[[1]][ii])/(predictionstest@tn[[1]][ii]+predictionstest@fp[[1]][ii]+predictionstest@fn[[1]][ii])
  FDRo[ii]<-predictionstest@fp[[1]][ii]/(predictionstest@fp[[1]][ii]+predictionstest@tp[[1]][ii])
}
R1 <- rocdemo.sca(rbinom(40,1,.3), rnorm(40), dxrule.sca, 
                  caseLabel="new case", markerLabel="demo Marker" )
R1@sens<-TPRo; #Sensitivity
R1@spec<-(1-FPRo);#Specificity 
ACCo<-(TPRo+TNRo)/(TPRo+TNRo+FPRo+FNRo)
MCCo<-((TPRo*TNRo)-(FPRo*FNRo))/sqrt((TPRo+FPRo)*(TPRo+FNRo)*(TNRo+FPRo)*(TNRo+FNRo))
MCRo<-(FPRo+FNRo)/(TPRo+TNRo+FPRo+FNRo)
FDRo<-(FPRo)/(TPRo+FPRo)
FNRo=(FNRo)/(FNRo+TPRo)
FDRo = (FPRo/(TPRo+FPRo))
AUCo<-print(AUC(R1));   # AUC
pAUCo<-print(pAUC(R1, 0.1))   # pAUC upto FPR<=0.1
FPRo
TPRo[[16]]
TNRo[[16]]
FNRo[[16]]
FPRo[[16]]
FDRo[[16]]
ACCo[[16]]
MCCo[[16]]
MCRo[[16]]


##
###############
tiff("ROC_RF_All_PC_Train_Test_Final.tiff", width = 5.5, height = 4.5, res = 300, units = "in", compression = c("lzw"))
roc.curve(data_test1$class, Pred.RFF_test_1[,1], col = "red", lwd = '1.5', main="ROC curve for RF Classifier in Glioblastoma")
roc.curve(data_test2$class, Pred.RFF_test_2[,1], col = "blue", lwd = '1.5',add = TRUE)
legend(.22, .22, legend = c("AUC:0.992(Test Data)", "AUC:0.989 (Independent test Data)"),
       col = c("red","blue"), lty = 1, lwd = 1.5)
dev.off()


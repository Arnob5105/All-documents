##==========================================================##
##************* Random Forest ******************************##
#===========================================================##
rm(list = ls())
getwd()
setwd("F:\\0.Paper_1_PC\\Paper1_Old\\9. Classification\\Group_1") 
#install.packages("randomForest")
#install.packages("ROCR")
#install.packages("ROSE")
#install.packages("ggplot2")
#install.packages("pROC")
require(randomForest)
require(ROCR)
require(ROSE)
require(ggplot2)
require(pROC)

data_train<- read.csv(file="train_1.csv", header = TRUE, sep=",")
data_test<- data.frame(read.csv(file = "test_1.csv", header = TRUE, sep = ","))
Data_Normal_H  <- read.csv(file="Data_Normal_H.csv", header = TRUE, sep=",")
Hub <- Data_Normal_H[,1]
#View(data_train)
class(data_train$X)

########## Train the data ###########
#install.packages("e1071")
#install.packages("caret")
#install.packages("rpart")
#install.packages("gtools")
require("e1071")
require("rpart")
require("caret")
require("gtools")
data_train$X = as.factor(data_train$X)
class(data_train$X)


#set.seed(100)
rfFit <-  randomForest(data_train$X ~ ADAM10+COL1A1+COL1A2+COL3A1+FBN1+FN1
                       +LAMC1+P4HB ,data = data_train, 
                        ntree = 100, norm.votes=FALSE)
rfFit
#plot(rfFit)

#### Training data set #########
Pred.RFF_tr <- predict(rfFit, data_train, "prob")
head(Pred.RFF_tr)
table(predict(rfFit),data_train$X)
predictionstr = prediction(Pred.RFF_tr[,2], data_train$X) # Tumor
predictionstr
performance(predictionstr,"tpr","fpr")



####Confusion matrix####
m=length(predictionstr@cutoffs[[1]]);m
pDEo<-rep(0,m); TPRo=rep(0,m); TNRo=rep(0,m); FPRo=rep(0,m); FNRo=rep(0,m); 
ERo=rep(0,m); FDRo=rep(0,m);ACCRo=rep(0,m);MCCRo=rep(0,m);MCRRo=rep(0,m);  
for (ii in 1:m) 
{
  TPRo[ii]<-predictionstr@tp[[1]][ii]/(predictionstr@tp[[1]][ii]+predictionstr@fn[[1]][ii])
  FPRo[ii]<-predictionstr@fp[[1]][ii]/(predictionstr@tn[[1]][ii]+predictionstr@fp[[1]][ii])
  TNRo[ii]<-predictionstr@tn[[1]][ii]/(predictionstr@tn[[1]][ii]+predictionstr@fp[[1]][ii])
  FNRo[ii]<-predictionstr@fn[[1]][ii]/(predictionstr@fn[[1]][ii]+predictionstr@tp[[1]][ii])
  ACCRo[ii]<-(predictionstr@tp[[1]][ii]+predictionstr@tn[[1]][ii])/(predictionstr@tn[[1]][ii]+predictionstr@fp[[1]][ii]+predictionstr@fp[[1]][ii]+predictionstr@fn[[1]][ii])
  MCCRo[ii]<-((predictionstr@tp[[1]][ii]*predictionstr@tn[[1]][ii])-(predictionstr@fp[[1]][ii]*predictionstr@fn[[1]][ii]))/sqrt((predictionstr@tp[[1]][ii]+predictionstr@fp[[1]][ii])*(predictionstr@tp[[1]][ii]+predictionstr@fn[[1]][ii])*(predictionstr@tn[[1]][ii]+predictionstr@fp[[1]][ii])*(predictionstr@tn[[1]][ii]+predictionstr@fn[[1]][ii]))
  MCRRo[ii]<-(predictionstr@fp[[1]][ii]+predictionstr@fn[[1]][ii])/(predictionstr@tn[[1]][ii]+predictionstr@fp[[1]][ii]+predictionstr@fn[[1]][ii])
  FDRo[ii]<-predictionstr@fp[[1]][ii]/(predictionstr@fp[[1]][ii]+predictionstr@tp[[1]][ii])
}

##=========================
#BioConductor installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(version = "3.12")
##=========================
#BiocManager::install("ROC")
citation("ROC")
require(ROC)
R1 <- rocdemo.sca(rbinom(40,1,.3), rnorm(40), dxrule.sca, 
                  caseLabel="new case", markerLabel="demo Marker" )
R1@sens<-TPRo; #Sensitivity (TPR)
R1@spec<-(1-FPRo);#Specificity (TNR)
ACC1<-(TPRo+TNRo)/(TPRo+TNRo+FPRo+FNRo) # Accuracy of model
MCC1<-((TPRo*TNRo)-(FPRo*FNRo))/sqrt((TPRo+FPRo)*(TPRo+FNRo)*(TNRo+FPRo)*(TNRo+FNRo))
MCR1<-(FPRo+FNRo)/(TPRo+TNRo+FPRo+FNRo)
FDR1<-(FPRo)/(TPRo+FPRo)
FNR1=(FNRo)/(FNRo+TPRo)
FDR = (FPRo/(TPRo+FPRo))
AUC1<-print(AUC(R1));   # AUC
pAUC1<-print(pAUC(R1, 0.2))   # pAUC upto FPR<=0.2
FPRo
TPRo[[30]]
TNRo[[30]]
FPRo[[30]]
FNRo[[30]]
FDR1[[30]]
ACC1[[30]]
MCC1[[30]]
MCR1[[30]]
AUC(R1)
### Confusion matrix
library(caret)
confusionMatrix(as.factor(data_train$X), predict(rfFit))


###TEST SET####
Pred.RFF_test <- predict(rfFit, data_test, "prob")
table(predict(rfFit, data_test), data_test$X)

predictionstest = prediction(Pred.RFF_test[,2], data_test$X)
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
ACC1<-(TPRo+TNRo)/(TPRo+TNRo+FPRo+FNRo)
MCC1<-((TPRo*TNRo)-(FPRo*FNRo))/sqrt((TPRo+FPRo)*(TPRo+FNRo)*(TNRo+FPRo)*(TNRo+FNRo))
MCR1<-(FPRo+FNRo)/(TPRo+TNRo+FPRo+FNRo)
FDR1<-(FPRo)/(TPRo+FPRo)
FNR1=(FNRo)/(FNRo+TPRo)
FDR = (FPRo/(TPRo+FPRo))
AUC1<-print(AUC(R1));   # AUC
pAUC1<-print(pAUC(R1, 0.2))   # pAUC upto FPR<=0.2
FPRo
TPRo[[17]]
TNRo[[17]]
FNRo[[17]]
FPRo[[17]]
FDR1[[17]]
ACC1[[17]]
MCC1[[17]]
MCR1[[17]]

###ROC####
tiff("ROC_RF_L.tiff", width = 5.5, height = 4.5, res = 300, units = "in", compression = c("lzw"))
roc.curve(data_train$X, Pred.RFF_tr[,1], col = "red", lwd = '1.1',
          main="ROC curve for RF Classifier")
text(.7,.3, "AUC : 1.00 (Train data)", col = 'red')
roc.curve(data_test$X, Pred.RFF_test[,1], col = "green", lwd = '1.1', add = TRUE)
text(.7,.2, "AUC : 0.945  (Test data)", col='green')
#segments(0, 0, 1, 1, lty=2)
#Axis(side=1, at=c(0.0, .20, .40, .60, .80, 1.0), labels = c("0%", "20%", "40%", "60%", "80%", "100%"))
#Axis(side=2, at=c(0.0, .20, .40, .60, .80, 1.0), labels = c("0%", "20%", "40%", "60%", "80%", "100%"))
dev.off()


##==============================================================
## ********************* SVM ***********************************
##==============================================================
##directory set
library(e1071)
library(ROCR)
library(ROSE)
#library(ROC)

data_train <- read.csv(file="train_1.csv", header = TRUE, sep=",")
data_test <- data.frame(read.csv(file = "test_1.csv", header = TRUE, sep = ","))

#### Train Model #########
status_s = as.factor(data_train$X)
svmmodel<- svm(status_s ~ ADAM10+COL1A1+COL1A2+COL3A1+FBN1+FN1
               +LAMC1+P4HB, 
               data = data_train)
svmmodel

###For Train Set########
svmmodel.predicttr<-predict(svmmodel,data_train ,decision.values=TRUE)
svmmodel.probstr<-attr(svmmodel.predicttr,"decision.values")
predictionstr<-prediction(svmmodel.probstr,data_train$X)

####Confusion matrix train ####
m=length(predictionstr@cutoffs[[1]]);m
pDEo<-rep(0,m); TPRo=rep(0,m); TNRo=rep(0,m); FPRo=rep(0,m); FNRo=rep(0,m); 
ERo=rep(0,m); FDRo=rep(0,m);ACCRo=rep(0,m);MCCRo=rep(0,m);MCRRo=rep(0,m);  
for (ii in 1:m) 
{
  TPRo[ii]<-predictionstr@tp[[1]][ii]/(predictionstr@tp[[1]][ii]+predictionstr@fn[[1]][ii])
  FPRo[ii]<-predictionstr@fp[[1]][ii]/(predictionstr@tn[[1]][ii]+predictionstr@fp[[1]][ii])
  TNRo[ii]<-predictionstr@tn[[1]][ii]/(predictionstr@tn[[1]][ii]+predictionstr@fp[[1]][ii])
  FNRo[ii]<-predictionstr@fn[[1]][ii]/(predictionstr@fn[[1]][ii]+predictionstr@tp[[1]][ii])
  ACCRo[ii]<-(predictionstr@tp[[1]][ii]+predictionstr@tn[[1]][ii])/(predictionstr@tn[[1]][ii]+predictionstr@fp[[1]][ii]+predictionstr@fp[[1]][ii]+predictionstr@fn[[1]][ii])
  MCCRo[ii]<-((predictionstr@tp[[1]][ii]*predictionstr@tn[[1]][ii])-(predictionstr@fp[[1]][ii]*predictionstr@fn[[1]][ii]))/sqrt((predictionstr@tp[[1]][ii]+predictionstr@fp[[1]][ii])*(predictionstr@tp[[1]][ii]+predictionstr@fn[[1]][ii])*(predictionstr@tn[[1]][ii]+predictionstr@fp[[1]][ii])*(predictionstr@tn[[1]][ii]+predictionstr@fn[[1]][ii]))
  MCRRo[ii]<-(predictionstr@fp[[1]][ii]+predictionstr@fn[[1]][ii])/(predictionstr@tn[[1]][ii]+predictionstr@fp[[1]][ii]+predictionstr@fn[[1]][ii])
}
R1 <- rocdemo.sca(rbinom(40,1,.3), rnorm(40), dxrule.sca, 
                  caseLabel="new case", markerLabel="demo Marker" )
R1@sens<-TPRo; #Sensitivity
R1@spec<-(1-FPRo);#Specificity 
R1@cuts<-predictionstr@cutoffs[[1]];
ACC<-(TPRo+TNRo)/(TPRo+TNRo+FPRo+FNRo)
MCC<-((TPRo*TNRo)-(FPRo*FNRo))/sqrt((TPRo+FPRo)*(TPRo+FNRo)*(TNRo+FPRo)*(TNRo+FNRo))
MCR<-(FPRo+FNRo)/(TPRo+TNRo+FPRo+FNRo)
FDR = (FPRo/(TPRo+FPRo))
AUC1<-print(AUC(R1)); # AUC
pAUC1<-print(pAUC(R1,0.2))# pAUC upto FPR<=0.2
FPRo
TPRo[[]]
TNRo[[]]
FNRo[[]]
ACC[[]]
MCC[[]]
MCR[[]]
FDR[[]]
####*******Confusion matrix**************
table(predict(svmmodel),data_train$X,dnn=c("Prediction","Actual"))
library(caret)
confusionMatrix(as.factor(data_train$X), predict(svmmodel))


###For Test Set########
svmmodel.predicttest<-predict(svmmodel,data_test ,decision.values=TRUE)
svmmodel.probstest<-attr(svmmodel.predicttest,"decision.values")
svmmodel.classtest<-predict(svmmodel,data_test, type="class")
predictionstest<-prediction(svmmodel.probstest,data_test$X)

#roc analysis for test data
predictionstest<-prediction(svmmodel.probstest,data_test$X)
svmmodel.performance_SVM <- performance(predictionstest,"tpr","fpr")
svmmodel.auc<-performance(predictionstest,"auc")@y.values[[1]];svmmodel.auc


####Confusion matrix test ####
m=length(predictionstest@cutoffs[[1]]);m
pDEo<-rep(0,m); TPRo=rep(0,m); TNRo=rep(0,m); FPRo=rep(0,m); FNRo=rep(0,m); 
ERo=rep(0,m); FDRo=rep(0,m);ACCRo=rep(0,m);MCCRo=rep(0,m);MCRRo=rep(0,m);  
for (ii in 1:m) 
{
  TPRo[ii]<-predictionstest@tp[[1]][ii]/(predictionstest@tp[[1]][ii]+predictionstest@fn[[1]][ii])
  FPRo[ii]<-predictionstest@fp[[1]][ii]/(predictionstest@tn[[1]][ii]+predictionstest@fp[[1]][ii])
  TNRo[ii]<-predictionstest@tn[[1]][ii]/(predictionstest@tn[[1]][ii]+predictionstest@fp[[1]][ii])
  FNRo[ii]<-predictionstest@fn[[1]][ii]/(predictionstest@fn[[1]][ii]+predictionstest@tp[[1]][ii])
  ACCRo[ii]<-(predictionstest@tp[[1]][ii]+predictionstest@tn[[1]][ii])/(predictionstest@tn[[1]][ii]+predictionstest@fp[[1]][ii]+predictionstest@fp[[1]][ii]+predictionstest@fn[[1]][ii])
  MCCRo[ii]<-((predictionstest@tp[[1]][ii]*predictionstest@tn[[1]][ii])-(predictionstest@fp[[1]][ii]*predictionstest@fn[[1]][ii]))/sqrt((predictionstest@tp[[1]][ii]+predictionstest@fp[[1]][ii])*(predictionstest@tp[[1]][ii]+predictionstest@fn[[1]][ii])*(predictionstest@tn[[1]][ii]+predictionstest@fp[[1]][ii])*(predictionstest@tn[[1]][ii]+predictionstest@fn[[1]][ii]))
  MCRRo[ii]<-(predictionstest@fp[[1]][ii]+predictionstest@fn[[1]][ii])/(predictionstest@tn[[1]][ii]+predictionstest@fp[[1]][ii]+predictionstest@fn[[1]][ii])
}
R1 <- rocdemo.sca(rbinom(40,1,.3), rnorm(40), dxrule.sca, 
                  caseLabel="new case", markerLabel="demo Marker" )
R1@sens<-TPRo; #Sensitivity
R1@spec<-(1-FPRo);#Specificity 
R1@cuts<-predictionstest@cutoffs[[1]];
ACC<-(TPRo+TNRo)/(TPRo+TNRo+FPRo+FNRo)
MCC<-((TPRo*TNRo)-(FPRo*FNRo))/sqrt((TPRo+FPRo)*(TPRo+FNRo)*(TNRo+FPRo)*(TNRo+FNRo))
MCR<-(FPRo+FNRo)/(TPRo+TNRo+FPRo+FNRo)
FDR = (FPRo/(TPRo+FPRo))
AUC1<-print(AUC(R1)); # AUC
pAUC1<-print(pAUC(R1,0.2))# pAUC upto FPR<=0.2
FPRo
TPRo[[6]]
TNRo[[4]]
FNRo[[4]]
ACC[[4]]
MCC[[10]]
MCR[[10]]
FDR[[10]]

###ROC####
tiff("ROC_SVM_L.tiff", width = 5.5, height = 4.5, res = 300, units = "in", compression = c("lzw"))
roc.curve(data_train$X, svmmodel.probstr[,1], col = "red", lwd = '1.1',
          main="ROC curve for SVM Classifier")
text(.7,.3, "AUC : 0.953   (Train data)", col = 'red')
roc.curve(data_test$X, svmmodel.probstest[,1], col = "green", lwd = '1.1', add = TRUE)
text(.7,.2, "AUC : 0.955  (Test data)", col='green')
#segments(0, 0, 1, 1, lty=2)
#Axis(side=1, at=c(0.0, .20, .40, .60, .80, 1.0), labels = c("0%", "20%", "40%", "60%", "80%", "100%"))
#Axis(side=2, at=c(0.0, .20, .40, .60, .80, 1.0), labels = c("0%", "20%", "40%", "60%", "80%", "100%"))
dev.off()

##==============================================================
## ********************** AdaBost ******************************
##==============================================================
library("fastAdaboost")
library(ROCR)
library(ROC)
library(ROSE)
library(pROC)
library(ggplot2)
library(plotROC)

data_train<- read.csv(file="train_1.csv", header=TRUE, sep=",")
data_test<- data.frame(read.csv(file="test_1.csv", header=TRUE, sep=","))

####For Training Set###
data_train$X =factor(data_train$X)
#model.adaboost<-adaboost(X ~., data_train , 20)
model.adaboost<-adaboost(X ~ADAM10+COL1A1+COL1A2+COL3A1+FBN1+FN1
                         +LAMC1+P4HB,data = data_train , 15)

######## For training datasets ###########
predADA <- predict( model.adaboost, newdata=data_train)
predictions_t = prediction(predADA$prob[,2], data_train$X)


###confusion Matrix####
m=length(predictions_t@cutoffs[[1]]);m
pDEo<-rep(0,m); TPRo=rep(0,m); TNRo=rep(0,m); FPRo=rep(0,m); FNRo=rep(0,m); 
ERo=rep(0,m); FDRo=rep(0,m);ACCRo=rep(0,m);MCCRo=rep(0,m);MCRRo=rep(0,m); 
for (ii in 1:m) 
{
  TPRo[ii]<-predictions_t@tp[[1]][ii]/(predictions_t@tp[[1]][ii]+predictions_t@fn[[1]][ii])
  FPRo[ii]<-predictions_t@fp[[1]][ii]/(predictions_t@tn[[1]][ii]+predictions_t@fp[[1]][ii])
  TNRo[ii]<-predictions_t@tn[[1]][ii]/(predictions_t@tn[[1]][ii]+predictions_t@fp[[1]][ii])
  FNRo[ii]<-predictions_t@fn[[1]][ii]/(predictions_t@fn[[1]][ii]+predictions_t@tp[[1]][ii])
  ACCRo[ii]<-(predictions_t@tp[[1]][ii]+predictions_t@tn[[1]][ii])/(predictions_t@tn[[1]][ii]+predictions_t@fp[[1]][ii]+predictions_t@fp[[1]][ii]+predictions_t@fn[[1]][ii])
  MCCRo[ii]<-((predictions_t@tp[[1]][ii]*predictions_t@tn[[1]][ii])-(predictions_t@fp[[1]][ii]*predictions_t@fn[[1]][ii]))/sqrt((predictions_t@tp[[1]][ii]+predictions_t@fp[[1]][ii])*(predictions_t@tp[[1]][ii]+predictions_t@fn[[1]][ii])*(predictions_t@tn[[1]][ii]+predictions_t@fp[[1]][ii])*(predictions_t@tn[[1]][ii]+predictions_t@fn[[1]][ii]))
  MCRRo[ii]<-(predictions_t@fp[[1]][ii]+predictions_t@fn[[1]][ii])/(predictions_t@tn[[1]][ii]+predictions_t@fp[[1]][ii]+predictions_t@fn[[1]][ii])
}
R1 <- rocdemo.sca(rbinom(40,1,.3), rnorm(40), dxrule.sca, 
                  caseLabel="new case", markerLabel="demo Marker" )
R1@sens<-TPRo; #Sensitivity
R1@spec<-(1-FPRo);#Specificity 
R1@cuts<-predictions_t@cutoffs[[1]];
ACC1<-(TPRo+TNRo)/(TPRo+TNRo+FPRo+FNRo)
MCC1<-((TPRo*TNRo)-(FPRo*FNRo))/sqrt((TPRo+FPRo)*(TPRo+FNRo)*(TNRo+FPRo)*(TNRo+FNRo))
MCR1<-(FPRo+FNRo)/(TPRo+TNRo+FPRo+FNRo)
FDR = (FPRo/(TPRo+FPRo))
AUC1<-print(AUC(R1)); # AUC
pAUC1<-print(pAUC(R1,0.2))# pAUC upto FPR<=0.2

FPRo
TPRo[[65]]
FNRo[[65]]
TNRo[[65]]  
ACC1[[65]]
MCC1[[65]]
MCR1[[65]]
FDR[[65]]

####For Test Set#####
pred_test = predict(model.adaboost, newdata=data_test)
predictionstest = prediction(pred_test$prob[,2], data_test$X)

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
  ACCRo[ii]<-(predictionstest@tp[[1]][ii]+predictionstest@tn[[1]][ii])/(predictionstest@tn[[1]][ii]+predictionstest@fp[[1]][ii]+predictionstest@fp[[1]][ii]+predictionstest@fn[[1]][ii])
  MCCRo[ii]<-((predictionstest@tp[[1]][ii]*predictionstest@tn[[1]][ii])-(predictionstest@fp[[1]][ii]*predictionstest@fn[[1]][ii]))/sqrt((predictionstest@tp[[1]][ii]+predictionstest@fp[[1]][ii])*(predictionstest@tp[[1]][ii]+predictionstest@fn[[1]][ii])*(predictionstest@tn[[1]][ii]+predictionstest@fp[[1]][ii])*(predictionstest@tn[[1]][ii]+predictionstest@fn[[1]][ii]))
  MCRRo[ii]<-(predictionstest@fp[[1]][ii]+predictionstest@fn[[1]][ii])/(predictionstest@tn[[1]][ii]+predictionstest@fp[[1]][ii]+predictionstest@fn[[1]][ii])
}
R1 <- rocdemo.sca(rbinom(40,1,.3), rnorm(40), dxrule.sca, 
                  caseLabel="new case", markerLabel="demo Marker" )
R1@sens<-TPRo; #Sensitivity
R1@spec<-(1-FPRo);#Specificity 
R1@cuts<-predictionstest@cutoffs[[1]];
ACC1<-(TPRo+TNRo)/(TPRo+TNRo+FPRo+FNRo)
MCC1<-((TPRo*TNRo)-(FPRo*FNRo))/sqrt((TPRo+FPRo)*(TPRo+FNRo)*(TNRo+FPRo)*(TNRo+FNRo))
MCR1<-(FPRo+FNRo)/(TPRo+TNRo+FPRo+FNRo)
FDR = (FPRo/(TPRo+FPRo))
AUC1<-print(AUC(R1)); # AUC
pAUC1<-print(pAUC(R1,0.2))# pAUC upto FPR<=0.2
FPRo
TPRo[[12]]
FNRo[[12]]
TNRo[[12]]
ACC1[[12]]
MCC1[[12]]
MCR1[[12]]
FDR[[12]]

###ROC####
tiff("ROC_Adaboost.tiff", width = 5.5, height = 4.5, res = 300, units = "in", compression = c("lzw"))

roc.curve(data_train$X, predADA$prob[,2],col = "red", lwd = "1.2",main="ROC curve for Adaboost")
text(.73, 0.35, "AUC : 1.00 (train data)", col = "red")

roc.curve(data_test$X,pred_test$prob[,2],col = 'green',add = TRUE, lwd = '1.2')
text(.73, 0.24, "AUC : 0.886 (test data)", col = "green")
dev.off()


##=========================================================
## ****************** Naive Bayes *************************
##=========================================================

library(e1071)
library(ROSE)
library(ROC)
library(ROCR)
library(pROC)

data_train<- read.csv(file="train_1.csv", header = TRUE, sep=",")
data_test<- data.frame(read.csv(file = "test_1.csv", header = TRUE, sep = ","))


set.seed(100)
model.NB <- naiveBayes(X ~ ADAM10+COL1A1+COL1A2+COL3A1+FBN1+FN1
                       +LAMC1+P4HB, data = data_train)
model.NB

### Training data #######
Pred.trbinNB <- predict(model.NB, data_train)
Pred.NB.trbin <- predict(model.NB, data_train,"raw")
forestpredtrbin = prediction(Pred.NB.trbin[,2], data_train$X)
svmmodel.performanc_NB <- performance(forestpredtrbin,"tpr","fpr")


### Confusion Matrix for training ####
m=length(forestpredtrbin@cutoffs[[1]]);m
pDEo<-rep(0,m); TPRo=rep(0,m); TNRo=rep(0,m); FPRo=rep(0,m); FNRo=rep(0,m); 
ERo=rep(0,m); FDRo=rep(0,m);ACCRo=rep(0,m);MCCRo=rep(0,m);MCRRo=rep(0,m);  
for (ii in 1:m) 
{
  TPRo[ii]<-forestpredtrbin@tp[[1]][ii]/(forestpredtrbin@tp[[1]][ii]+forestpredtrbin@fn[[1]][ii])
  FPRo[ii]<-forestpredtrbin@fp[[1]][ii]/(forestpredtrbin@tn[[1]][ii]+forestpredtrbin@fp[[1]][ii])
  TNRo[ii]<-forestpredtrbin@tn[[1]][ii]/(forestpredtrbin@tn[[1]][ii]+forestpredtrbin@fp[[1]][ii])
  FNRo[ii]<-forestpredtrbin@fn[[1]][ii]/(forestpredtrbin@fn[[1]][ii]+forestpredtrbin@tp[[1]][ii])
  ACCRo[ii]<-(forestpredtrbin@tp[[1]][ii]+forestpredtrbin@tn[[1]][ii])/(forestpredtrbin@tn[[1]][ii]+forestpredtrbin@fp[[1]][ii]+forestpredtrbin@fn[[1]][ii]+forestpredtrbin@tp[[1]][ii])
  MCCRo[ii]<-((forestpredtrbin@tp[[1]][ii]*forestpredtrbin@tn[[1]][ii])-(forestpredtrbin@fp[[1]][ii]*forestpredtrbin@fn[[1]][ii]))/(sqrt((forestpredtrbin@tp[[1]][ii]+forestpredtrbin@fp[[1]][ii])*(forestpredtrbin@tp[[1]][ii]+forestpredtrbin@fn[[1]][ii])*(forestpredtrbin@tn[[1]][ii]+forestpredtrbin@fp[[1]][ii])*(forestpredtrbin@tn[[1]][ii]+forestpredtrbin@fn[[1]][ii])))
  MCRRo[ii]<-(forestpredtrbin@fp[[1]][ii]+forestpredtrbin@fn[[1]][ii])/(forestpredtrbin@tn[[1]][ii]+forestpredtrbin@fp[[1]][ii]+forestpredtrbin@fn[[1]][ii]+forestpredtrbin@tp[[1]][ii])
}
R1 <- rocdemo.sca(rbinom(40,1,.3), rnorm(40), dxrule.sca, 
                  caseLabel="new case", markerLabel="demo Marker" )
R1@sens<-TPRo; #Sensitivity
R1@spec<-(1-FPRo);#Specificity 
R1@cuts<-forestpredtrbin@cutoffs[[1]];
Accuracy<-(TPRo+TNRo)/(TPRo+TNRo+FPRo+FNRo)
MCC<-((TPRo*TNRo)-(FPRo*FNRo))/sqrt((TPRo+FPRo)*(TPRo+FNRo)*(TNRo+FPRo)*(TNRo+FNRo))
MCR<-(FPRo+FNRo)/(TPRo+TNRo+FPRo+FNRo)
FDR = (FPRo/(TPRo+FPRo))
AUC1<-print(AUC(R1)); # AUC
pAUC1<-print(pAUC(R1,0.2))# pAUC upto FPR<=0.2
FPRo
TPRo[[14]]
TNRo[[14]]
FNRo[[14]]
Accuracy[[14]]
MCC[[14]]
MCR[[14]]
FDR[[14]]

### Test data #######
Pred.testbinNB <- predict(model.NB, data_test)
Pred.NB.testbin <- predict(model.NB, data_test,"raw")
forestpredtestbin = prediction(Pred.NB.testbin[,2], data_test$X)
svmmodel.performanc_NB <- performance(forestpredtestbin,"tpr","fpr")

###confusion Matrix####
m=length(forestpredtestbin@cutoffs[[1]]);m
pDEo<-rep(0,m); TPRo=rep(0,m); TNRo=rep(0,m); FPRo=rep(0,m); FNRo=rep(0,m); 
ERo=rep(0,m); FDRo=rep(0,m);ACCRo=rep(0,m);MCCRo=rep(0,m);MCRRo=rep(0,m);  
for (ii in 1:m) 
{
  TPRo[ii]<-forestpredtestbin@tp[[1]][ii]/(forestpredtestbin@tp[[1]][ii]+forestpredtestbin@fn[[1]][ii])
  FPRo[ii]<-forestpredtestbin@fp[[1]][ii]/(forestpredtestbin@tn[[1]][ii]+forestpredtestbin@fp[[1]][ii])
  TNRo[ii]<-forestpredtestbin@tn[[1]][ii]/(forestpredtestbin@tn[[1]][ii]+forestpredtestbin@fp[[1]][ii])
  FNRo[ii]<-forestpredtestbin@fn[[1]][ii]/(forestpredtestbin@fn[[1]][ii]+forestpredtestbin@tp[[1]][ii])
  ACCRo[ii]<-(forestpredtestbin@tp[[1]][ii]+forestpredtestbin@tn[[1]][ii])/(forestpredtestbin@tn[[1]][ii]+forestpredtestbin@fp[[1]][ii]+forestpredtestbin@fn[[1]][ii]+forestpredtestbin@tp[[1]][ii])
  MCCRo[ii]<-((forestpredtestbin@tp[[1]][ii]*forestpredtestbin@tn[[1]][ii])-(forestpredtestbin@fp[[1]][ii]*forestpredtestbin@fn[[1]][ii]))/(sqrt((forestpredtestbin@tp[[1]][ii]+forestpredtestbin@fp[[1]][ii])*(forestpredtestbin@tp[[1]][ii]+forestpredtestbin@fn[[1]][ii])*(forestpredtestbin@tn[[1]][ii]+forestpredtestbin@fp[[1]][ii])*(forestpredtestbin@tn[[1]][ii]+forestpredtestbin@fn[[1]][ii])))
  MCRRo[ii]<-(forestpredtestbin@fp[[1]][ii]+forestpredtestbin@fn[[1]][ii])/(forestpredtestbin@tn[[1]][ii]+forestpredtestbin@fp[[1]][ii]+forestpredtestbin@fn[[1]][ii]+forestpredtestbin@tp[[1]][ii])
}
R1 <- rocdemo.sca(rbinom(40,1,.3), rnorm(40), dxrule.sca, 
                  caseLabel="new case", markerLabel="demo Marker" )
R1@sens<-TPRo; #Sensitivity
R1@spec<-(1-FPRo);#Specificity 
R1@cuts<-forestpredtestbin@cutoffs[[1]];
Accuracy<-(TPRo+TNRo)/(TPRo+TNRo+FPRo+FNRo)
MCC<-((TPRo*TNRo)-(FPRo*FNRo))/sqrt((TPRo+FPRo)*(TPRo+FNRo)*(TNRo+FPRo)*(TNRo+FNRo))
MCR<-(FPRo+FNRo)/(TPRo+TNRo+FPRo+FNRo)
FDR = (FPRo/(TPRo+FPRo))
AUC1<-print(AUC(R1)); # AUC
pAUC1<-print(pAUC(R1,0.2))# pAUC upto FPR<=0.2
FPRo
TPRo[[29]]
TNRo[[29]]
FNRo[[29]]
Accuracy[[29]]
MCC[[29]]
MCR[[29]]
FDR[[29]]

##ROC###
tiff("ROC curve for NB classifier.tiff", width = 5.5, height = 4.5, res = 300, units = "in", compression = c("lzw"))
roc.curve(data_train$X,Pred.NB.trbin[,1],col = 'red', lwd = '1.5',
          main="ROC curve for Naive Bayes Classifier")
text(.7,.3, "AUC : 0.833 (train data)", col = "red")
roc.curve(data_test$X,Pred.NB.testbin[,1],col = 'green', lwd = '1.5', add = TRUE)
text(.7,.2, "AUC : 0.915 (test data)", col = "green")
dev.off()




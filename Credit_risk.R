# LOAD THE DATA

# R code by Tyler Moore for analyzing Bitcoin exchange failure
# Read in survival table
survex3<-read.table("https://tylermoore.ens.utulsa.edu/data/bitcoin/sdatamodCC.csv", sep=",",header=T)
survex3$dailyvol<-survex3$totalvol/survex3$lifetime
# Anti-money laundering indicators
aml<-read.table("https://tylermoore.ens.utulsa.edu/data/bitcoin/compliance-aml-cft-whole.csv",
                sep=",",header=T)
# Merge the two dataset by country
amls<-merge(survex3,aml,by="Country")
# Prepare the variables for the survival analysis
amls$PrevN<-amls$Preventive/amls$PreventativeMax
amls$InstN<-amls$Institutional/amls$InstitutionalMax
amls$LegalN<-amls$Legal/amls$LegalMax
amls$AllN<-amls$All/amls$AllMax
amlsv<-amls[!is.na(amls$dailyvol),]
row.names(amlsv)<-amlsv$exchange

# express operational life in months
amlsv$lifemo<-amlsv$lifetime/30


# extract only the data we need from the previous amlsv dataset
cleandata<-amlsv[, c("censored", "dailyvol", "lifemo", "Hacked", "All")]

# Logit model for the PD of bitcoin exchanges
bit.logit<- glm(censored ~ dailyvol + lifemo + Hacked + All, data = cleandata,
                family = binomial(link="logit"))
summary(bit.logit)

## Call: ...
## Deviance Residuals:
##     Min       1Q   Median       3Q      Max
## -1.5735  -0.8528  -0.4719   0.8691   1.9753
## Coefficients:
##               Estimate Std. Error z value Pr(>|z|)
## (Intercept)  3.411e+00  2.125e+00   1.605   0.1085
## dailyvol    -4.014e-05  1.291e-04  -0.311   0.7558
## lifemo      -1.346e-01  5.722e-02  -2.352   0.0187 *
## HackedTRUE   8.941e-01  9.633e-01   0.928   0.3533
## All         -7.770e-02  6.948e-02  -1.118   0.2635
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## (Dispersion parameter for binomial family taken to be 1)
##     Null deviance: 55.051  on 39  degrees of freedom
## Residual deviance: 44.768  on 35  degrees of freedom
## AIC: 54.768
## Number of Fisher Scoring iterations: 5



bit.logit2<-glm(censored~lifemo,data=cleandata,family=binomial(link="logit"))
summary(bit.logit2)

## Call: ...
## glm(formula = censored ~ lifemo, family = binomial(link = "logit"),
##     data = cleandata)
## Deviance Residuals:
##     Min       1Q   Median       3Q      Max
## -1.7926  -0.8446  -0.5697   0.8456   1.9370
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)
## (Intercept)   1.4245     0.6993   2.037  0.04165 *
## lifemo       -0.1389     0.0532  -2.610  0.00904 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## (Dispersion parameter for binomial family taken to be 1)
##     Null deviance: 55.051  on 39  degrees of freedom
## Residual deviance: 46.605  on 38  degrees of freedom
## AIC: 50.605


range(cleandata$lifemo)
#[1]  0.30000 30.46667

#create a sequence of values for lifemo to produce fitted values
xlife <- seq(0, 31, 0.1)
ylife<-predict(bit.logit2,data.frame(lifemo=xlife),type="response",se.fit=TRUE)
plot(cleandata$lifemo, cleandata$censored, pch = 16,
 xlab = "EXCHANGE OPERATIONAL LIFE IN MONTHS", ylab = "PROBABILITY OF DEFAULT")
lines(xlife, ylife$fit)
# Add standard errors
lines(xlife, ylife$fit - (1.96*ylife$se.fit), lty=2)
lines(xlife, ylife$fit + (1.96*ylife$se.fit), lty=2)

library(margins)
margins::margins(bit.logit2)
#Average marginal effects
#glm(formula=censored~lifemo,family=binomial(link="logit"),data=cleandata)
#   lifemo
# -0.02756
cplot(bit.logit2, "lifemo")


library(MASS)
bit.lda<-MASS::lda(censored~dailyvol+lifemo+Hacked+All,data=cleandata);bit.lda
## Call: ...
## Prior probabilities of groups:
##    0    1
## 0.55 0.45
## Group means:
##      dailyvol    lifemo HackedTRUE      All
##   0 2593.5731 15.121212  0.1818182 28.84164
##   1  642.7426  8.468519  0.2777778 27.47783
## Coefficients of linear discriminants:
##                      LD1
##   dailyvol   -1.260748e-05
##   lifemo     -1.327506e-01
##   HackedTRUE  7.972208e-01
##   All        -7.186039e-02


plot(bit.lda, col = as.integer(cleandata$censored))


# Decision tree for the PD of bitcoin exchanges
library(rpart)
bit.DT<- rpart(censored ~ dailyvol + lifemo + Hacked + All, data = cleandata,
                method="class")
bit.DT
# n= 40
# node), split, n, loss, yval, (yprob),       * denotes terminal node
#   1) root 40 18 0 (0.5500000 0.4500000)
#   2) lifemo>=5.266667 30  9 0 (0.7000000 0.3000000) *
#   3) lifemo< 5.266667 10  1 1 (0.1000000 0.9000000) *

# Let's compare the basic plot (left), with a more informative plot (right)
par(mfrow=c(1,2))
library(rattle); library(rpart.plot); library(RColorBrewer);
plot(bit.DT); text(bit.DT);    fancyRpartPlot(bit.DT)

# Random Forest for the PD of bitcoin exchanges
library(randomForest)
# The response variables should be a factor to have a classification,
# otherwise a regression is assumed
bit.RF<- randomForest(as.factor(censored) ~ dailyvol + lifemo + Hacked + All,
      data = cleandata, importance = TRUE);  bit.RF
## Call: ...
##             Type of random forest: classification
##                      Number of trees: 500
## No. of variables tried at each split: 2
##         OOB estimate of  error rate: 42.5%
## Confusion matrix: # We will come later to this concept!
##    0 1 class.error
## 0 14 8   0.3636364
## 1  9 9   0.5000000

importance(bit.RF)
##                  0          1 MeanDecreaseAccuracy MeanDecreaseGini
## dailyvol -2.576670 -3.7100668            -4.064215        5.4260828
## lifemo    9.544349  9.0971138            11.427514        9.3317095
## Hacked   -4.995400 -0.4946245            -4.222217        0.6545848
## All      -2.173418 -1.1937564            -2.250970        3.7718467

varImpPlot(bit.RF)


# Conditional inference Random Forest for the PD of bitcoin exchanges
library(party)
bit.CRF<- cforest(as.factor(censored) ~ dailyvol + lifemo + Hacked + All,
      data = cleandata, control = cforest_unbiased(mtry = 2))

# Conditional variable importance for ‘cforest’, following the permutation
# principle of the ‘mean decrease in accuracy’ importance in ‘randomForest’.
varimp(bit.CRF, conditional = TRUE)
#     dailyvol        lifemo        Hacked           All
# 0.0051428571  0.0717142857 -0.0005714286 -0.0088571429

# We can see an individual tree (from the ensemble of 500 trees) as follows
party:::prettytree(bit.CRF@ensemble[[3]],names(bit.CRF@data@get("input")))
## 
## 1) lifemo <= 6.333333; criterion = 0.992, statistic = 6.95
##   2)*  weights = 0
## 1) lifemo > 6.333333
##   3)*  weights = 0


# Support VEctor Machine for the PD of bitcoin exchanges
library(e1071)
bit.lsvm<-svm(censored ~ dailyvol + lifemo + Hacked + All,
      data = cleandata, type = 'C-classification',  kernel='linear')
bit.lsvm

## Call: ...
## Parameters:
##    SVM-Type:  C-classification
##  SVM-Kernel:  linear
##        cost:  1
##       gamma:  0.2
## Number of Support Vectors:  29

# % correctly predicted in the training sample with a linear SVM
predicted.lsvm <- predict(bit.lsvm, data= cleandata)
100*mean(predicted.lsvm==cleandata$censored)
# [1] 67.5
# % correctly predicted in the training sample with a nonlinear SVM
bit.nlsvm<-svm(censored ~ dailyvol + lifemo + Hacked + All,
      data = cleandata, type = 'C-classification', kernel='radial')
predicted.nlsvm <- predict(bit.nlsvm, data= cleandata)
100*mean(predicted.nlsvm==cleandata$censored)
## [1] 75

# Finally,we can try to improve the performance by doing some parameter tuning:
tune_par <- tune.svm(censored ~ dailyvol + lifemo + Hacked + All,data = cleandata,
                  gamma=10^(-3:3),cost=c(0.01,0.1,1,10,100,1000),kernel="radial")

bit.nlsvm2<-svm(censored ~ dailyvol + lifemo + Hacked + All,
  data=cleandata, type='C-classification',cost=tune_par$best.parameters$cost,
  gamma=tune_par$best.parameters$gamma, kernel='radial')

predicted.nlsvm2 <- predict(bit.nlsvm2, data= cleandata)
100*mean(predicted.nlsvm2==cleandata$censored)
#[1] 100

## # ==> of course, this is in-sample fitting! We'll see later how all these models
## # perform with an out-of-sample dataset.


## =============== MODEL EVALUATION ==================================

# 1) We first consider a DATA SPLIT 67%/33% with our clean dataset
library(caret); library(ROCR); split=0.67; set.seed(1)
trainIndex <- caret::createDataPartition(cleandata$censored, p=split, 
                                         list=FALSE)
data_train <- cleandata[ trainIndex,]
data_test <- cleandata[-trainIndex,]
# Train a logit model
train.logit<- glm(censored ~ dailyvol + lifemo + Hacked + All, 
                  data=data_train, family=binomial(link="logit"))
# Make forecasts
pred.logit <- predict(train.logit, data_test, type="response")

# 1A) Computing ROC plot and AUC using the ROCR package
pred <- ROCR::prediction(pred.logit , data_test$censored)
perf1A <- ROCR::performance(pred, measure="tpr", x.measure="fpr")
plot(perf1A)


auc <- ROCR::performance(pred, measure = "auc")
auc@y.values[[1]]


# 1B) Computing ROC plot and AUC using the pROC package
perf1B<- pROC::roc(censored ~ pred.logit, data = data_test)
perf1B$auc
plot(perf1B)


# If I want to find the cut-off with the maximum accuracy,then I do as follows:
pROC::coords(perf1B, "best")
pred.logit.best=ifelse(pred.logit>=0.7988907,1,0)
caret::confusionMatrix(as.factor(pred.logit.best), as.factor(data_test$censored))


# 2) We now consider a k-fold Cross Validation with our clean dataset
#    Specify training control  
set.seed(1)
train_control.k10 <- caret::trainControl(method="cv", number=10)
# train the logit model
train.logit.k10 <- caret::train(as.factor(censored) ~ dailyvol+lifemo+Hacked+All, 
                                data = data_train,trControl=train_control.k10, method="glm", 
                                family = binomial(link="logit"))
# summarize results
print(train.logit.k10)   

# 3) We now consider a REPEATED k-fold Cross Validation with our clean dataset
#    Specify training control  
set.seed(1)
train_control.Rk10 <- trainControl(method="repeatedcv", number=10, repeats=10)
# Train the logit model
train.logit.Rk10 <- train(as.factor(censored) ~ dailyvol + lifemo + Hacked + All,
                          data = data_train,trControl=train_control.Rk10, method="glm", 
                          family = binomial(link="logit"))
# Summarize results
print(train.logit.Rk10)   

# 4) Finally, we consider the Leave One Out Cross Validation (LOOCV)
#    Specify training control  
train_control.LOOCV <- trainControl(method="LOOCV")
# train the logit model
train.logit.LOOCV <- train(as.factor(censored) ~ dailyvol + lifemo + Hacked + All,
                           data = data_train, trControl=train_control.LOOCV, method="glm",
                           family = binomial(link="logit"))
# Summarize results
print(train.logit.LOOCV)   

# Optimal thresholds
library(OptimalCutpoints)
data.fore.logit<-data.frame(pred.logit=pred.logit,actual.data=data_test$censored)
opt.cut<-OptimalCutpoints::optimal.cutpoints(X="pred.logit",status="actual.data",
                                             tag.healthy = 0,methods = "MCT", data = data.fore.logit, 
                                             control = control.cutpoints(CFP=1, CFN=45))
summary(opt.cut)

opt.cut2<-OptimalCutpoints::optimal.cutpoints(X ="pred.logit",status ="actual.data",
                                             tag.healthy = 0,methods = "MCT", data = data.fore.logit,
                                             control = control.cutpoints(CFP=1, CFN=1))
summary(opt.cut2)


## comparison of the three machine learning models
# Prepare training scheme
control <- trainControl(method="repeatedcv", number=10, repeats=20)

# Train the nonlinear SVM model
set.seed(1)
modelsvmRadial <- train(as.factor(censored) ~ dailyvol + lifemo + Hacked + All,
                        data = cleandata,  method="svmRadial", trControl=control)
# Train the Random Forest model
set.seed(1)
modelrf <- train(as.factor(censored) ~ dailyvol + lifemo + Hacked + All, 
                 data = cleandata,  method="rf", trControl=control)
# Train the Random forest with conditional inference trees
set.seed(1)
modelcforest <- train(as.factor(censored) ~ dailyvol + lifemo + Hacked + All,
                      data = cleandata,  method="cforest", trControl=control)
# Collect resamples
results <- resamples(list(SVM=modelsvmRadial, RF=modelrf, CFOREST=modelcforest))
# Summarize the results
summary(results)


# Boxplots of results
bwplot(results)
# Dot-plots of results
dotplot(results)

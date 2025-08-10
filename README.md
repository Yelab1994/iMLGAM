# iMLGAM: integrated Machine Learning and Genetic Algorithm-driven Multiomics analysis for pan-Cancer immunotherapy response prediction

## Overview

This repository contains data and code related to the study of predictive models for pan-cancer immunotherapy. The purpose of this study is to construct a signature based on RNA sequencing data to predict the outcomes of immunotherapy in patients with pan-cancer.

## Directory Structure

- `Figure1`, `Figure2`, ...: Data and original code for saving generated figures.
- The source code of the iMLAGA package.

## Prerequisites

Make sure you have R and the following packages installed: iMLAGA. Install "iMLAGA" via this [GitHub page](https://github.com/Yelab1994/iMLAGA) or by running the code.

```R
devtools::install_github("Yelab1994/iMLGAM")
```
## Step 1: Constructing Base Learners
Download sample files via this [Baidu Netdisk page](https://pan.baidu.com/s/1hlTEJXgVYTQfrrMBXHe-RA?pwd=s8hr) and load the necessary R packages.

```R
#load the necessary R packages
library(qs)
library(abess)
library(iMLGAM)
#Read the example file
dat <- qread("Sample File.qs")
#The first column is the outcome (with row names as 'res'), and the subsequent columns are features.
for (i in names(dat)) {
  dat[[i]] <- dat[[i]][,-c(1,3,4)]
  colnames( dat[[i]])[1] <- "res"
}
#Determining the training set
trainset <- dat$trainset_7
#Determining the validation set
validationset <- dat$trainset_3
#Determining the test set 
testset <- dat$testset 
#Feature Selection using ABESS Algorithm
trainset <- trainset[!is.na(trainset$res),]
abess_fit <- abess(trainset[,-1], trainset$res,support.size=0:100)
key_gene <- extract(abess_fit, 5)[4][1]$support.vars
trainset <- trainset[,c("res",key_gene )]
validationset  <- validationset [,c("res",key_gene )]
testset <- testset [,c("res",key_gene )]
#Constructing Base Learners
basic_learner_list <- basic_learner(trainingset=trainset , CVnumber=10, Cvrepeats=5, ncore=8)
```

## Step 2: Constructing prediction matrices for various machine learning models
Batch run the prediction results of models constructed in Step 1, and merge them into a single matrix.
```R
#Obtaining the prediction results of each model in the training set
trainres <- basic_learner_predictor(trainset, basic_learner_list )
#Obtaining the prediction results of each model in the validation set
validatres <- basic_learner_predictor(validationset, basic_learner_list )
#Obtaining the prediction results of each model in the test set
testres <- basic_learner_predictor(testset, basic_learner_list )
```

## Step 3: Running a genetic learning algorithm to obtain the optimal combination of base learners
```R
#Merging with the prediction matrix
dat <- qread("ICI cohort gene pairs.qs")
for (i in names(dat)) {
  dat[[i]] <- dat[[i]][,-c(1,3,4)]
  colnames( dat[[i]])[1] <- "res"
}
jj <- intersect(rownames(dat$trainset_7),rownames(trainres))
trainres <- cbind(res=dat$trainset_7[jj,]$res,trainres[jj,])
jj <- intersect(rownames(dat$trainset_3),rownames(validatres))
validatres <- cbind(res=dat$trainset_3[jj,]$res,validatres[jj,])
jj <- intersect(rownames(dat$testset),rownames(testres))
testres <- cbind(res=dat$testset[jj,]$res,testres[jj,])
#Merge with the prediction matrix
ga_fit <- GA_train(trainres)
```

## Step 4: Calculating Risk Scores for Each set
Scoring each set using the formula constructed in Step 3.
```R
#Obtaining the evaluation score of the training set
trainres $riskscore <- predict(ga_fit ,trainres)
#Obtaining the evaluation score of the validation set
validatres $riskscore <- predict(ga_fit ,validatres)
#Obtaining the evaluation score of the test set
testres$riskscore<- predict(ga_fit ,testres)
```

setwd("D:\\科研\\遗传算法R包项目\\iMLAGA\\iMLGAM\\1.Figure1")
library(qs)
library(abess)
library(iMLGAM)

dat <- qread("ICI cohort gene pairs.qs")
for (i in names(dat)) {
  dat[[i]] <- dat[[i]][,-c(1,3,4)]
  colnames( dat[[i]])[1] <- "res"
}
trainset <- dat$trainset_7
validationset <- dat$trainset_3
testset <- dat$testset
###ABESS
trainset <- trainset[!is.na(trainset$res),]
abess_fit <- abess(trainset[,-1], trainset$res,support.size=0:100)
key_gene <- extract(abess_fit, 5)[4][1]$support.vars
trainset <- trainset[,c("res",key_gene )]
validationset  <- validationset [,c("res",key_gene )]
testset <- testset [,c("res",key_gene )]
#############step1 basic learn
basic_learner_list <- basic_learner(trainingset=trainset , CVnumber=10, Cvrepeats=5, ncore=64)
######step2 basic learn predict
trainres <- basic_learner_predictor(trainset, basic_learner_list )
validatres <- basic_learner_predictor(validationset, basic_learner_list )
trainres <- basic_learner_predictor(trainset, basic_learner_list )
testres <- basic_learner_predictor(testset, basic_learner_list )
####setp3 GA
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
ga_fit <- GA_train(trainres)
###riskscore
trainres $riskscore <- predict(ga_fit ,trainres)
validatres $riskscore <- predict(ga_fit ,validatres)
testres$riskscore<- predict(ga_fit ,testres)
roc(trainres$res,trainres $riskscore )
roc(testres$res,testres$riskscore )
roc(validatres$res,validatres$riskscore )
riskscore <- list()
riskscore$trainset_7 <- trainres
riskscore$trainset_3 <- validatres
riskscore$testset  <- testres
###output
saveRDS(riskscore,"riskscore.rds")
saveRDS(basic_learner_predictor,"basic_learner_predictor.rds")
saveRDS(ga_fit,"ga_fit.rds")
saveRDS(key_gene,"key_gene_pairs.rds")


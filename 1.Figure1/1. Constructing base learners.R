library(data.table)
library(tibble)
library(survival)
library(randomForestSRC)
library(survivalsvm)
library(gbm)
library(caret)
library(abess)
library(glmnet)
library(e1071)
library(kknn)
library(ggpubr)
library(ggthemes)
library(qs)
library(iMLGAM)
library(doParallel)
dat <- qread("ICI_data.qs")
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
#################step1 Constructing Base Learners
basic_learner_list <- basic_learner(trainingset=trainset , CVnumber=10, Cvrepeats=5, ncore=8)
################step2 Constructing prediction matrices for various machine learning models
#Obtaining the prediction results of each model in the training set
trainres <- basic_learner_predictor(trainset, basic_learner_list )
#Obtaining the prediction results of each model in the validation set
validatres <- basic_learner_predictor(validationset, basic_learner_list )
#Obtaining the prediction results of each model in the test set
testres <- basic_learner_predictor(testset, basic_learner_list )
################step3 Running the genetic algorith
trainres <- cbind(res=na.omit(dat$trainset_7$res),trainres)
validatres <- cbind(res=na.omit(dat$trainset_3$res),validatres)
testres  <- cbind(res=na.omit(dat$testset$res),testres )
ga_fit <- GA_train(trainres)
################step 4 Calculating Risk Scores for each set
#Obtaining the evaluation score of the training set
trainres $riskscore <- predict(ga_fit ,trainres)
#Obtaining the evaluation score of the validation set
validatres $riskscore <- predict(ga_fit ,validatres)
#Obtaining the evaluation score of the test set
testres$riskscore<- predict(ga_fit ,testres)
roc(trainres$res,trainres$riskscore)
roc(testres$res,testres$riskscore)
roc(validatres$res,validatres$riskscore)
#################step 5 visualization
###ROC
bas_learn <- list(trainset=trainres,testset=testres,validationset=validatres)
roc_list <- list()
boxplot_list <- list()
for (i in names(bas_learn)) {   dat <- bas_learn[[i]]
dat <- dat[!is.na(dat$res),]
roc_obj <- roc(dat$res, dat$riskscore,smooth=F) 
roc_coords <- coords(roc_obj)
roc_list [[i]]<- ggplot(data = roc_coords, aes(x = 1-specificity, y = sensitivity)) +
  geom_ribbon(aes(ymin = 0, ymax =  sensitivity), fill = "#029978", alpha = 0.2)+
  geom_line(color = "#029978",lwd=1)  +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black", size = 1)+
  labs(title = "ROC Curve",
       x = "1 - Specificity",
       y = "Sensitivity") +
  annotate("text",x=0.7,y=0.3,label=paste0("AUC = ",round(as.numeric(roc_obj $auc),3)),
           colour = "black", size = 5, fontface = 2)+
  theme(
    axis.text = element_text(face = "plain", colour = "black", size = 12, hjust = 0.5),
    axis.line = element_blank(),
    plot.title = element_text(face = "bold", colour = "black", size = 16, hjust = 0.5),
    title = element_text(face = "bold", colour = "black", size = 13),
    panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid"),
    panel.background = element_rect(fill = "#F1F6FC"),
    panel.grid.major = element_line(color = "#CFD3D6", size = .5, linetype = "dotdash"),
    legend.position = "none"
  )
#Boxplot
dat$res <- ifelse(dat$res==1,"Non-res","Res")
boxplot_list[[i]]  <- ggplot(data = dat, aes(x = res, y = riskscore, fill = res)) +     ##需要修改
  scale_fill_manual(values = c("#C31820", "#106C61")) +
  geom_boxplot(notch = FALSE, outlier.size = 1, color = "black", lwd = 0.8, alpha = 1,width=0.8,outlier.colour = NA) +
  theme_classic() +
  ylab('riskScore') +
  xlab("") +
  ggtitle(i) +
  theme(
    axis.text = element_text(face = "plain", colour = "black", size = 12, hjust = 0.5),
    axis.line = element_blank(),
    plot.title = element_text(face = "bold", colour = "black", size = 16, hjust = 0.5),
    title = element_text(face = "bold", colour = "black", size = 13),
    panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid"),
    panel.background = element_rect(fill = "#F1F6FC"),
    panel.grid.major = element_line(color = "#CFD3D6", size = .5, linetype = "dotdash"),
    legend.position = "none"
  ) +
  stat_compare_means(label.x =1.5,aes(group=res),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

}
####################################save
qsave(basic_learner_list ,"basic_learner_list.qs")
qsave(ga_fit ,"ga_fit.qs")

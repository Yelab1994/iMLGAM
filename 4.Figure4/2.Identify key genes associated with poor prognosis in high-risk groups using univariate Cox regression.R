setwd("D:\\科研\\ICBP_code\\Figure4")
library(data.table)
library(tibble)
library(qs)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(survival)
library(dplyr)
library(ggplot2)
###Read the output risk scores.
rt <- readRDS("riskscore.rds")
rt <- rbind(rt$trainset_7,rt$trainset_3,rt$testset)
rt $Risk <- ifelse(rt$riskscore>0.996,"High","Low") ##0.996 is the optimal cutoff value determined in our study
rt_id <- row.names(rt)[rt$Risk=="High"]
###Read ICI chort RNA seq
dat <- qread("ICI RNA-seq.qs")
dat <- rbind(dat$trainset_7,dat$trainset_3,dat$testset)
rt <- dat [rt_id,]
rt <- rt[,3:ncol(rt)]


###Identify key genes associated with poor prognosis in high-risk groups using univariate Cox regression
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
  dat <- rt[,c("futime","fustat",i)]
  colnames(dat)[3] <- "gene"
  cox <- coxph(Surv(futime, fustat) ~., data = dat)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients["gene","Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int["gene","exp(coef)"],
                     HR.95L=coxSummary$conf.int["gene","lower .95"],
                     HR.95H=coxSummary$conf.int["gene","upper .95"],
                     pvalue=coxSummary$coefficients["gene","Pr(>|z|)"])
  )
}
outTab$pvalue <- as.numeric(outTab$pvalue)


library(IMvigor210CoreBiologies)
library(IOBR)
library(survival)
##Read the phenotype file of the Braun cohort.
cli1 <- read.csv("cli.csv",header = T,row.names = 1)
cli <- fread("RCC-Braun_2020.Response-clinical.tsv",header = T,data.table = F)
rownames(cli) <- cli$patient_name
###Merge with the output risk scores.
rt <- readRDS("riskscore.rds")
rt <- rbind(rt$trainset_7,rt$trainset_3,rt$testset)
jj <- intersect(rownames(phenoData),rownames(rt))
phenoData <- phenoData[jj,]
rt <- rt[jj,]
rt <- cbind(res=rt$res,riskscore =rt$riskscore,phenoData)
####Plot ROC Curve
ROC1 <- roc(rt$res,rt$riskscore)
ROC2 <- roc(rt$res,rt$TMB)
palette <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF")  
line_types <- c(1, 2, 3) 
line_width <- 2  
roc(as.factor(rt$res), rt$riskscore, plot=TRUE, print.thres=F, print.auc=F, asp=NA, col=palette[1], lwd=line_width, lty=line_types[1])

roc(as.factor(rt$res), rt$TMB, plot=TRUE, print.thres=F, print.auc=F, asp=NA, add=TRUE, col=palette[2], lwd=line_width, lty=line_types[2])
legend('bottomright',
       c(paste0('AUC for riskscore: ', sprintf("%.03f", as.numeric(ROC1$auc))),
         paste0('AUC for TMB: ', sprintf("%.03f", as.numeric(ROC2$auc)))),
       col=palette, lwd=line_width, lty=line_types, bty='n')
###################Draw a multifactorial regression forest plot.
table(rt$ImmunoPhenotype)
rt1 <- rt[,c(1,2,11,12,13,22,124)]
colnames(rt1)[6] <- "Metastasis"
rt1$Metastasis <- ifelse(rt1$Metastasis=="PRIMARY","No","Yes")
mod <- glm(res ~ ., data = rt1, family = binomial)
forest_model(mod)

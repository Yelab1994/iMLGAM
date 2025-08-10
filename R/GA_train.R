GA_train <- function(trainingset) {
  library(parallel)
  library(dplyr)
  library(caret)
  library(doParallel)
  library(recipes)
  library(glmnet)
  library(MASS)
  library(pROC)
  library(snow)
  library(genalg)
  fitnessFunction <- function (bitstring) {
    if (sum(bitstring ) == 0){return(1)}  else {
      selectedColumns <- which(bitstring == 1)+1
      datset <-  trainingset [,c(1, selectedColumns)]
      fin_fit<-glm(res ~ ., data = datset, family = "binomial")
      fin_fit<-stepAIC(fin_fit, direction = "both")
      datset$riskscore <- predict( fin_fit ,datset)
      roc_value <- roc(datset$res,datset$riskscore)
      roc_value <- 1 - as.numeric(roc_value$auc)
      return(roc_value)
    }
  }
  monitor <- function(obj) {
    minEval = min(obj$evaluations);
    plot(obj, type="hist");
  }

  result <- rbga.bin(size=ncol(trainingset)-1, mutationChance=0.05, zeroToOneRatio=10,
                     evalFunc=fitnessFunction , verbose=F, monitorFunc=monitor,popSize=50, iters=20)

  combination <- result$population[which.min(result$evaluations), ]
  combination <- which(combination == 1)+1
  trainingset <- trainingset[,c(1,combination)]
  fin_fit<-glm(res ~ ., data = trainingset, family = "binomial")
  fin_fit<-stepAIC(fin_fit, direction = "both")
  return(fin_fit)
}


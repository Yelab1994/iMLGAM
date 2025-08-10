
basic_learner_predictor <- function(set,learnerlist){
  library(parallel)
  library(dplyr)
  library(caret)
  library(doParallel)
  library(recipes)
  library(glmnet)
  library(MASS)
  library(pROC)
  library(snow)
  ####predict ern
  ern_score <-  list()
  ern_fit <- learnerlist$ern_fit
  for (i in names(ern_fit)) {
    ern_score[[i]] <- predict(ern_fit[[i]]$fit, newdata = model.matrix( res~ ., set)[,-1], lambda = ern_fit[[i]]$lambda, type = "prob")[, 2]
  }
  ern_score <- do.call(cbind,ern_score)
  rownames(ern_score) <- rownames(model.matrix( res~ ., set)[,-1])
  ###predict KNN
  KNN_score <-  list()
  KNN_fit <- learnerlist$KNN_fit
  for (i in names(KNN_fit)) {
    KNN_score[[i]] <- predict(KNN_fit[[i]], newdata = model.matrix( res~ ., set)[,-1],  type = "prob")[, 2]
  }
  KNN_score <- do.call(cbind,KNN_score)
  rownames(KNN_score) <- rownames(model.matrix( res~ ., set)[,-1])
  ####predict RF
  rf_score <-  list()
  rf_fit <- learnerlist$rf_fit
  for (i in names(rf_fit)) {
    rf_score[[i]] <- predict(rf_fit[[i]], newdata = model.matrix( res~ ., set)[,-1],  type = "prob")[, 2]
  }
  rf_score <- do.call(cbind,rf_score)
  rownames(rf_score) <- rownames(model.matrix( res~ ., set)[,-1])
  ####predict svmLinear
  svmLinear_score <-  list()
  svmLinear_fit <- learnerlist$svmLinear_fit
  for (i in names(svmLinear_fit)) {
    svmLinear_score[[i]] <- as.numeric(predict(svmLinear_fit[[i]], newdata = model.matrix( res~ ., set)[,-1],type = "raw"))-1
  }
  svmLinear_score <- do.call(cbind,svmLinear_score)
  rownames(svmLinear_score) <- rownames(model.matrix( res~ ., set)[,-1])
  ###predict svmPoly
  svmPoly_score <-  list()
  svmPoly_fit <- learnerlist$svmPoly_fit
  for (i in names(svmPoly_fit)) {
    svmPoly_score[[i]] <- as.numeric(predict(svmPoly_fit[[i]], newdata = model.matrix( res~ ., set)[,-1],type = "raw"))-1
  }
  svmPoly_score <- do.call(cbind,svmPoly_score)
  rownames(svmPoly_score) <- rownames(model.matrix( res~ ., set)[,-1])
  ###predict svmRadial
  svmRadial_score <-  list()
  svmRadial_fit <- learnerlist$svmRadial_fit
  for (i in names(svmRadial_fit)) {
    svmRadial_score[[i]] <- as.numeric(predict(svmRadial_fit[[i]]$fit, newdata = model.matrix( res~ ., set)[,-1],type = "raw"))-1
  }
  svmRadial_score <- do.call(cbind,svmRadial_score)
  rownames(svmRadial_score) <- rownames(model.matrix( res~ ., set)[,-1])
  all_riskscore <- list(ern_score=ern_score,KNN_score=KNN_score,rf_score=rf_score,svmLinear_score=svmLinear_score,svmPoly_score=svmPoly_score,svmRadial_score=svmRadial_score)
  all_riskscore <-as.data.frame( do.call(cbind,all_riskscore ))
  return(all_riskscore)
}

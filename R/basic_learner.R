
basic_learner <- function(trainingset,CVnumber,Cvrepeats,ncore) {
  library(parallel)
  library(dplyr)
  library(caret)
  library(doParallel)
  library(recipes)
  library(glmnet)
  library(MASS)
  library(pROC)
  library(snow)
  cl <- makePSOCKcluster(ncore)
  registerDoParallel(cl)
  colnames(trainingset)[1] <- "res"
  #trainControl parameters
  train_control <- trainControl(method = "repeatedcv",
                                number = CVnumber,
                                repeats = Cvrepeats,
                                verboseIter = FALSE,
                                sampling = "smote")
  ####ERN
  alpha <- c(0,0.2,0.4,0.6,0.8,1.0)
  x <- model.matrix( res~ ., trainingset)[, -1]
  y <- as.factor(trainingset$res)
  ern_fit <- list()
  for (i in alpha  ) {

    model <- train(
      x, y,
      method = "glmnet",
      tuneGrid = expand.grid(alpha = i, lambda=seq(0.001, 0.1, length=100)),
      trControl = train_control
    )
    lambda <- model[["bestTune"]][["lambda"]]
    id <- paste0("ERN_",which(i==alpha))
    ern_fit[[id]]$"fit"<- model
    ern_fit[[id]]$"lambda" <- lambda
  }
  ###RF
  depth <- c(3,5,10 ,20)
  ntree <- c(200, 500, 800 ,1200)
  rf_fit <- list()
  for (i in depth ) {for(yy in ntree ){
    model <- train(
      x, y,
      method = "rf",
      trControl = train_control,
      ntree=yy,
      depth=i
    )
    id <- paste0("RSF_",which(i==depth)*4-(4-which(yy==ntree)))
    rf_fit[[id]] <- model
  }
  }
  ####svmLinear
  parm = c(1, 10, 100, 1000)
  kernel = c("svmLinear", "svmRadial", "svmPoly")

  ####svmRadial
  svmRadial_fit <- list()
  for (i in parm  ) {
    model <- train(
      x, y,
      method = "svmRadial",
      trControl = train_control,
      tuneGrid = expand.grid(sigma = seq(0.001, 0.1, length = 10),
                             C =i)
    )
    id <- paste0("svmRadial_",which(i==parm ))
    sigma=model[["bestTune"]][["sigma"]]

    svmRadial_fit [[id]]$"fit" <- model
    svmRadial_fit[[id]]$"sigma" <- sigma

  }

  ####svmLinear
  svmLinear_fit <- list()
  for (i in parm  ) {
    model <- train(
      x, y,
      method = "svmLinear",
      trControl = train_control,
      tuneGrid = expand.grid(
        C=i)
    )

    id <- paste0("svmLinear_",which(i==parm ))
    svmLinear_fit [[id]] <- model

  }
  ####svmPoly
  svmPoly_fit <- list()
  for (i in parm  ) {
    model <- train(
      x, y,
      method = "svmPoly",
      trControl = train_control,
      tuneGrid = expand.grid(degree=3, scale=1,
                             C=i)
    )

    id <- paste0("svmPoly_",which(i==parm ))
    svmPoly_fit[[id]] <- model

  }
  ################KNN
  k = c(3, 5, 7, 9)
  kernel = c("optimal", "rank", "rectangular", "triangular")
  KNN_fit <- list()
  for (i in k  ) { for(yy in kernel){
    model <- train(
      x, y,
      method = "kknn",
      trControl = train_control,
      tuneGrid = expand.grid(kmax=i,distance = c(1),kernel =yy)
    )
    id <- paste0("KNN_",which(i==k)*4-(4-which(yy==kernel)))
    KNN_fit [[id]] <- model

  }
  }
  basic_learner_list <- list(ern_fit=ern_fit,KNN_fit=KNN_fit,
                             rf_fit=rf_fit,svmLinear_fit=svmLinear_fit,
                             svmPoly_fit=svmPoly_fit,svmRadial_fit=svmRadial_fit)

  return(basic_learner_list)

}

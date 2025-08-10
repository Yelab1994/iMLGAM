# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

gene_pair_matrix <- function(dat,combination,ncore){
  library(parallel)
  library(dplyr)
  library(caret)
  library(doParallel)
  library(recipes)
  library(glmnet)
  library(MASS)
  library(pROC)
  library(snow)
  cl <- makeCluster(ncore)
  clusterExport(cl, as.character(substitute(dat)))
  clusterExport(cl, as.character(substitute(combination)))
  clusterExport(cl,varlist = c("gene_pair"),envir = environment())
  res <- parLapply(cl, 1:ncol(combination), function(i) {
    gene_pair(dat, combination[1,i], combination[2,i])
  })
  stopCluster(cl)
  res <- do.call(cbind,res)
  return(res)
}


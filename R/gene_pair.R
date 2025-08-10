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

gene_pair <- function(dat, gene1, gene2) {
  library(parallel)
  library(dplyr)
  library(caret)
  library(doParallel)
  library(recipes)
  library(glmnet)
  library(MASS)
  library(pROC)
  library(snow)
  dat_combination <- dat[, c(gene1, gene2)]
  pair <- ifelse(as.numeric(dat_combination[, 1]) > as.numeric(dat_combination[, 2]), 1, 0)
  pair <- as.matrix(pair)
  colnames(pair) <- paste0(colnames(dat_combination)[1], "|", colnames(dat_combination)[2])
  rownames(pair) <- rownames(dat_combination)
  return(pair)
}



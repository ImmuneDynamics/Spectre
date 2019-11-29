library(class)
library(caret)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
source('../R/classify.R')

data <- demo.clustered
classifier.col <- c(5,6,8,9,11,13,17:19,21:29,32)
label.col <- c(39)

knn_model <- classify(data, classifier.col, label.col, num.neighbours = 5)

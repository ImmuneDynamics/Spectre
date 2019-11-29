library(class)
library(Spectre)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('../R/classify.R')

data <- demo.clustered
classifier.col <- c(5,6,8,9,11,13,17:19,21:29,32)
label.col <- c(39)

# fine tune your number of neighbours
accuracy <- train.classifier(data, classifier.col, label.col, num.neighbours = 5)

unlabelled.data <- demo.start[1:180, classifier.col]
training.data <- data[, classifier.col]
training.label <- data[,label.col]

predicted.label <- run.classifier(training.data = training.data,
               training.label = training.label,
               unlabelled.data = unlabelled.data,
               num.neighbours = 5)

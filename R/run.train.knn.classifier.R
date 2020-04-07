#' train.knn.classifier - ...
#'
#' @usage train.knn.classifier(x, ...)
#'
#' @param train.data NO DEFAULT. A dataframe containing cells (rows) vs features/markers (columns) to be used to train classifier. 
#' Only include columns/features/markers to be used to train the classifier i.e. the ones that are useful in identifying cells (which population it belongs to).
#' @param label NO DEFAULT. A character vector containing the label (population name) for each cell in the train.data. Must be "CHARACTER" vector.
#' @param min.num.neighbours DEFAULTS to 1. When using a k-nearest neighbour (KNN) classifier, this parameter specifies the minimum number of nearest neighbours used to train the KNN classifier.
#' @param max.num.neighbours DEFAULTS to 1. When using a k-nearest neighbour (KNN) classifier, this parameter specifies the maximum number of nearest neighbours used to train the KNN classifier.
#' @param evaluation.metric DEFAULTS to 'Accuracy'. How do you want the classifier performance to be evaluated? By measuring accuracy (pass "Accuracy") or Cohen Kappa score (pass "Kappa").
#' @param num.folds DEFAULTS to 10. Number of chunks the training data is to be split into. The classifier will then take 1 chunk for testing the classifier (after it's trained), and use the remaining chunk to train the classifier.
#' @param num.repeats DEFAULTS to 10. Number of time the classifier will be trained (per number of neighbours). For each repeat, different chunk will be used for testing and training.
#' 
#' Train a k-nearest neighbour (KNN) classifier.
#' The classifier will be trained on a number of neighbours, starting from min.num.neighbours, and increased gradually by 1 until max.num.neighbours is reached.
#' For each number of neighbour, the accuracy of the classifier will be computed.
#' Data will be normalised to range between 0 and 1 before used in training the classifier.
#' The normalisation method used is MinMaxScaling method.
#' NOTE: the larger the training data, the longer it takes to train the classifier. Hence please be mindful of the training data size.
#'
#' @return The performance of the classifier on different number of neighbours as well as some description on how the training was performed
#' i.e. what sampling process is used (how the data is split for training and testing),
#' the sample sizes (number of data points used for training/testing).
#' Included as well is the recommended number of neighbours for the data and the reasoning behind why that number of neighbours is the best.
#'
#' @export

train.knn.classifier <- function(train.data, 
                                 label, 
                                 min.num.neighbours = 1,
                                 max.num.neighbours = 10,
                                 evaluation.metric = 'Accuracy',
                                 num.folds = 10,
                                 num.repeats = 10){
  
  # normalised the training data so each column is in range 0 to 1
  preprocess.model <- caret::preProcess(train.data, method = "range")
  normalised.data <- as.data.table(predict(preprocess.model, train.data))
  
  # add the label as a column in normalised data
  normalised.data[, ('population') := label]
  
  # split data into 10 folds and repeat training and testing 10 times
  # each iteration will train and test on different fold.
  trControl <- trainControl(
    method = 'repeatedcv',
    number = num.folds,
    repeats = num.repeats
  )
  
  # train the classifier
  knn.fit <- train(population ~ .,
                   method     = "knn",
                   tuneGrid   = expand.grid(k = min.num.neighbours:max.num.neighbours),
                   trControl  = trControl,
                   metric     = evaluation.metric,
                   data       = normalised.data)
  
  return(knn.fit)
}

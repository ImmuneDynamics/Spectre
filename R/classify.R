#' train.classifier - ...
#'
#' @usage train.classifier(x, ...)
#'
#' @param data NO DEFAULT. A dataframe containing cells (rows) vs features/markers (columns) that is the data to be classified.
#' @param classifier.col NO DEFAULT to KNN. Name of the column in the training data to be used as the training labels
#' @param label.col NO DEFAULT. Name of the new column that contains classification labels.
#' @param num.neighbours DEFAULTS to 1. When using a k-nearest neighbour classifier, then this parameter specifies the number of nearest neighbours.
#' @param training NO DEFAULT. A dataframe containing cells (rows) vs features/markers (columns) that is the data to be use as training.
#' 
#' Train a classifier and return the accuracy of the classifier.
#'
#' @return The accuracy of the classifier
#'
#' @export

train.classifier <- function(data, 
                     classifier.col,
                     label.col,
                     num.neighbours = 1){
  
  # split data into training and testing. TODO: convet to 10 fold cross validation
  ran <- sample(1:nrow(data), 0.9 * nrow(data)) 
  
  train.data <- data[ran, classifier.col]
  train.label <- data[ran, label.col]
  
  test.data <- data[-ran, classifier.col]
  test.label <- data[-ran, label.col]
  
  # min max normalisation function
  nor <-function(x) { (x -min(x))/(max(x)-min(x))   }
  
  train.data.norm <- as.data.frame(lapply(train.data[,], nor))
  test.data.norm <- as.data.frame(lapply(test.data[,], nor))
  
  pr <- knn(train.data.norm, test.data.norm, cl=train.label,k=num.neighbours)
  
  tab <- table(pr,test.label)
  
  ##this function divides the correct predictions by total number of predictions that tell us how accurate teh model is.
  
  accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}
  
  return(accuracy(tab))
}

#' run.classifier - ...
#'
#' @usage run.classifier(x, ...)
#'
#' @param training.data NO DEFAULT. A dataframe containing cells (rows) vs features/markers (columns) that is the data to be used to train classifer.
#' @param training.label NO DEFAULT to KNN. Name of the column in the training data to be used as the training labels.
#' @param unlabelled.data NO DEFAULT. A dataframe containing cells (rows) vs features/markers (columns) that is the data classified.
#' @param num.neighbours DEFAULTS to 1. When using a k-nearest neighbour classifier, then this parameter specifies the number of nearest neighbours.
#' 
#' Train and run the classifier.
#'
#' @return The predicted label for the unlabelled data.
#'
#' @export
run.classifier <- function(training.data, 
                           training.label,
                           unlabelled.data,
                           num.neighbours = 1){
  
  # min max normalisation function
  nor <-function(x) { (x -min(x))/(max(x)-min(x))   }
  
  # normalised training and unlabelled data together
  all.data <- rbind(training.data, unlabelled.data)
  all.data.norm <- as.data.frame(lapply(all.data, nor))
  
  train.data.norm <- all.data.norm[c(1:nrow(training.data)),]
  start.test.data <- nrow(training.data) + 1
  
  test.data.norm <- all.data.norm[c(start.test.data:nrow(all.data.norm)),]
  
  pr <- knn(train.data.norm, test.data.norm, cl=training.label,k=num.neighbours)
  
  return(pr)
}


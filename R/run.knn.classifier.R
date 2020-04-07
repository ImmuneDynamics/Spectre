#' run.knn.classifier - ...
#'
#' @usage run.knn.classifier(x, ...)
#'
#' @param train.data NO DEFAULT. A dataframe containing cells (rows) vs features/markers (columns) to be used to train classifier. 
#' Only include columns/features/markers to be used to train the classifier i.e. the ones that are useful in identifying cells (which population it belongs to).
#' @param train.label NO DEFAULT. A character vector containing the label (population name) for each cell in the train.data. Must be "CHARACTER" vector.
#' @param unlabelled.data NO DEFAULT. A dataframe containing cells (rows) vs features/markers (columns) that is the data classified.
#' @param num.neighbours DEFAULTS to 1. When using a k-nearest neighbour classifier, then this parameter specifies the number of nearest neighbours.
#' 
#' Train a k-NN classifier on a training data, and use it to classify an unlabelled data. 
#' Note that for the classifier to work as intended, the unlabelled.data has to be normalised to the range of the train.data.
#' Note make sure that train.data and unlabelled.data has exactly the same features/markers.
#'
#' @return The predicted label for the unlabelled data.
#'
#' @export
run.knn.classifier <- function(train.data, 
                               train.label,
                               unlabelled.data,
                               num.neighbours = 1){
  
  # normalised the training and unlabelled data so each column is in range 0 to 1
  # the unlabelled data will be normalised using the model built on the training data
  preprocess.model <- caret::preProcess(train.data, method = "range")
  norm.train.data <- as.data.table(predict(preprocess.model, train.data))
  norm.unlab.data <- as.data.table(predict(preprocess.model, unlabelled.data))
  
  pr <- knn(norm.train.data, norm.unlab.data, cl=train.label,k=num.neighbours)
  
  return(pr)
}

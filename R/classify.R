#' classify
#'
#' @usage classify(x, ...)
#'
#' @param x NO DEFAULT. A dataframe containing cells (rows) vs features/markers (columns) that is the data to be classified.
#' @param type DEFAULTS to KNN. Specify the type of classifier to be used. Can be KNN (k-nearest neighbour).
#' @param num.neighbours DEFAULTS to 1. When using a k-nearest neighbour classifier, then this parameter specifies the number of nearest neighbours.
#' @param training NO DEFAULT. A dataframe containing cells (rows) vs features/markers (columns) that is the data to be use as training.
#' @param training.label NO DEFAULT. Name of the column in the training data to be used as the training labels
#' @param new.label NO DEFAULT. Name of the new column that contains classification labels.
#'
#' @return ...
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @examples classify()
#'
#' @export

# MY OWN NOTE:
# Let's start with KNN for now.
# We need a dataframe for training data (we will split it into 10 fold stratified cross validation)
# Need the dataframe to classify.

classify <- function(data, 
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
  
  print(accuracy(tab))
  
  return(pr)
}


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

classify <- function(x,
                     type,
                     num.neighbours = 1,
                     training,
                     training.label,
                     new.label){
}

# MY OWN NOTE:
# Need method to just train and validate classifier
train.classifier <- function(training.data, training.label, num.neighbours) {
  
}


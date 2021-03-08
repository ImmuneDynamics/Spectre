#' run.knn.classifier - Run KNN classifier
#'
#' @usage run.knn.classifier(train.dat, unlabelled.dat, use.cols, label.col, num.neighbours)
#'
#' @param train.dat NO DEFAULT. data.frame. A dataframe containing cells (rows) vs features/markers (columns) to be used to train classifier. 
#' @param unlabelled.dat NO DEFAULT. data.frame. A dataframe containing cells (rows) vs features/markers (columns) that is the data classified.
#' @param use.cols NO DEFAULT. Vector of column names to use for training k-nearest neighbour (KNN) classifier.
#' @param label.col NO DEFAULT. Character. Name of the column representing the population name of the cells in dat.
#' @param num.neighbours DEFAULTS to 1. Numeric. When using a k-nearest neighbour classifier, then this parameter specifies the number of nearest neighbours.
#' 
#' Train a k-NN classifier on a training data, and use it to classify an unlabelled data. 
#' Note that for the classifier to work as intended, the unlabelled.data has to be normalised to the range of the train.data.
#' Note make sure that train.data and unlabelled.data has exactly the same features/markers.
#'
#' @return The predicted label for the unlabelled data.
#' 
#'
#'
#' @author Givanna Putri, \email{ghar1821@@uni.sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @export

run.knn.classifier <- function(train.dat,
                               unlabelled.dat,
                               use.cols,
                               label.col,
                               num.neighbours = 1,
                               seed=42){
  
  if(!is.element('caret', installed.packages()[,1])) stop('caret is required but not installed')
  if(!is.element('FNN', installed.packages()[,1])) stop('FNN is required but not installed')
  
  require(caret)
  require(FNN)
  
  set.seed(seed)
  
  # for testing
  # library(Spectre)
  # library(data.table)
  # unlabelled.dat <- do.subsample(demo.clustered, rep(10,6), seed = 20, divide.by = 'Population')[,c(11:19,25)]
  # train.dat <- demo.clustered[,c(11:19,25)]
  # use.cols <- names(demo.clustered)[c(11:19)]
  # # remove unlabelled dat that is in train.dat
  # train.dat <- train.dat[!unlabelled.dat, on=use.cols]
  # label.col <- 'Population'
  
  # isolate the features
  train.dat.features <- train.dat[, ..use.cols]
  unlabelled.dat.features <- unlabelled.dat[, ..use.cols]
  
  # isolate the label
  train.dat.labels <- as.character(train.dat[[label.col]])
  
  # normalised the training and unlabelled data so each column is in range 0 to 1
  # the unlabelled data will be normalised using the model built on the training data
  preprocess.model <- caret::preProcess(train.dat.features, method = "range")
  norm.train.data <- data.table(predict(preprocess.model, train.dat.features))
  norm.unlab.data <- data.table(predict(preprocess.model, unlabelled.dat.features))
  
  norm.train.data_mat <- as.matrix(norm.train.data)
  norm.unlab.data_mat <- as.matrix(norm.unlab.data)
  
  pr <- FNN::knn(train=norm.train.data_mat, 
                 test=norm.unlab.data_mat, 
                 cl=train.dat.labels,
                 k=num.neighbours)
  
  # append the predicted class
  predicted.cell.dat <- data.table(unlabelled.dat)
  predicted.cell.dat$Prediction <- pr
  
  # this gives the index of the nearest k neighbours data points
  closest_neighbours_indices <- attr(pr, "nn.index")
  closest_neighbours_indices <- data.table(closest_neighbours_indices)
  names(closest_neighbours_indices) <- sapply(c(1:ncol(closest_neighbours_indices)), function(x) paste0("Neighbour_", x))
  
  predicted.cell.dat <- cbind(predicted.cell.dat, closest_neighbours_indices)
  
  return(predicted.cell.dat)
}

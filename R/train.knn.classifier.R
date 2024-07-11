#' train.knn.classifier - Train KNN classifier
#'
#' @usage train.knn.classifier(dat, use.cols, label.col,min.num.neighbours, 
#' max.num.neighbours, evaluation.metric, num.folds, num.repeats)
#'
#' @param dat NO DEFAULT. data.frame. A dataframe containing cells (rows) vs features/markers (columns) to be used to train k-nearest neighbour (KNN) classifier
#' @param use.cols NO DEFAULT. Vector of column names to use for training k-nearest neighbour (KNN) classifier.
#' @param label.col NO DEFAULT. Character. Name of the column representing the population name of the cells in dat.
#' @param min.num.neighbours DEFAULTS to 1. Numeric. When using a k-nearest neighbour (KNN) classifier, this parameter specifies the minimum number of nearest neighbours used to train the KNN classifier.
#' @param max.num.neighbours DEFAULTS to 1. Numeric. When using a k-nearest neighbour (KNN) classifier, this parameter specifies the maximum number of nearest neighbours used to train the KNN classifier.
#' @param method DEFAULTS to random. Can either be random which randomly shuffle data and split them into 2 halves,
#'   one for training, the other for testing. Or CV (cross validation) where data is split into num.folds complementary
#'   portions, and num.folds-1 portions used for training and the remaining for testing.
#' @param num.folds DEFAULTS to 10. Numeric. Number of chunks the training data is to be split into. The classifier will then take 1 chunk for testing the classifier (after it's trained), and use the remaining chunk to train the classifier.
#' @param num.repeats DEFAULTS to 1. Numeric. Number of time the training data will be split into num.folds chunks.
#'   If you set this to 3 and num.folds to 10, the classifier will split data into 10 chunks, train classifier on those chunks 10 times,
#'   and repeat the entire procedure 3 times (each time differet data will be in each chunk).
#' @param seed DEFAULTS to 42. Seed used when splitting data into training and testing set.
#'
#' Train a k-nearest neighbour (KNN) classifier.
#' The classifier will be trained on a number of neighbours, starting from min.num.neighbours, and increased gradually by 1 until max.num.neighbours is reached.
#' For each number of neighbour, the accuracy (or other evaluation metric as specified in evaluation.metric) of the classifier will be computed.
#' Data will be normalised to range between 0 and 1 before used in training the classifier.
#' The normalisation method used is MinMaxScaling method.
#' NOTE: the larger the training data, the longer it takes to train the classifier. Hence please be mindful of the training data size.
#'
#' @return The performance of the classifier on different number of neighbours as well as some description on how the training was performed
#' i.e. what sampling process is used (how the data is split for training and testing),
#' the sample sizes (number of data points used for training/testing).
#' Included as well is the recommended number of neighbours for the data and the reasoning behind why that number of neighbours is the best.
#'
#'
#' @author Givanna Putri
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @export train.knn.classifier

train.knn.classifier <- function(dat,
                                 use.cols,
                                 label.col,
                                 min.num.neighbours = 1,
                                 max.num.neighbours = 10,
                                 method = c("random", "CV"),
                                 num.folds = 10,
                                 num.repeats = 1,
                                 seed = 42) {
  if (!is.element("caret", installed.packages()[, 1])) stop("caret is required but not installed")
  if (!is.element("data.table", installed.packages()[, 1])) stop("data.table is required but not installed")

  require(data.table)
  require(caret)

  dat.features <- dat[, ..use.cols]
  dat.labels <- as.character(dat[[label.col]])

  # normalised the training data so each column is in range 0 to 1
  preprocess.model <- caret::preProcess(dat.features, method = "range")
  normalised.data <- as.data.table(predict(preprocess.model, dat.features))

  # add the label as a column in normalised data
  normalised.data[, ("population") := dat.labels]

  method <- match.arg(method)

  set.seed(seed)

  if (method == "random") {
    num_row <- nrow(normalised.data)
    normalised.data <- normalised.data[sample(num_row), ]
    split <- round(num_row / 2)

    train.dat <- normalised.data[c(1:split), ]
    valid.dat <- normalised.data[c((split + 1):num_row), ]

    knn.res <- list()
    num_neighbours <- c(min.num.neighbours:max.num.neighbours)

    for (n in num_neighbours) {
      message(paste("Evaluating", n))
      res.knn <- Spectre::run.knn.classifier(
        train.dat = train.dat,
        unlabelled.dat = valid.dat,
        use.cols = use.cols,
        label.col = "population",
        num.neighbours = n
      )
      knn.res[[n]] <- res.knn
    }
    accuracies <- sapply(knn.res, function(res) {
      cm <- as.matrix(table(
        Actual = res$population,
        Predicted = res$Prediction
      ))
      corr.pred <- sum(diag(cm)) / length(res$Prediction)
    })

    acc.dat <- data.table(
      k = num_neighbours,
      accuracy = accuracies
    )
  } else if (method == "CV") {
    # split data into 10 folds and repeat training and testing 10 times
    # each iteration will train and test on different fold.
    trControl <- caret::trainControl(
      method = "repeatedcv",
      number = num.folds,
      repeats = num.repeats
    )

    # train the classifier
    knn.fit <- caret::train(population ~ .,
      method     = "knn",
      tuneGrid   = expand.grid(k = min.num.neighbours:max.num.neighbours),
      trControl  = trControl,
      metric     = "Accuracy",
      data       = normalised.data
    )

    acc.dat <- data.table(
      k = knn.fit$results$k,
      accuracy = knn.fit$results$Accuracy,
      accuracySD = knn.fit$results$AccuracySD
    )
  }
  return(acc.dat)
}

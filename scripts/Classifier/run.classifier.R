run.classifier <- function(train.dat,
                           unlabelled.dat,
                           use.cols,
                           label.col,
                           type = "knn",
                           knn.neighbours = 1){
  
  if(type == "knn"){
    
    if(!is.element('caret', installed.packages()[,1])) stop('caret is required but not installed')
    if(!is.element('class', installed.packages()[,1])) stop('class is required but not installed')
    
    require(caret)
    require(class)
    
    # isolate the features
    train.dat.features <- train.dat[, ..use.cols]
    unlabelled.dat.features <- unlabelled.dat[, ..use.cols]
    
    # isolate the label
    train.dat.labels <- as.character(train.dat[[label.col]])
    
    # normalised the training and unlabelled data so each column is in range 0 to 1
    # the unlabelled data will be normalised using the model built on the training data
    preprocess.model <- caret::preProcess(train.dat.features, method = "range")
    norm.train.data <- as.data.table(predict(preprocess.model, train.dat.features))
    norm.unlab.data <- as.data.table(predict(preprocess.model, unlabelled.dat.features))
    
    pr <- class::knn(norm.train.data, norm.unlab.data, cl=train.dat.labels,k=knn.neighbours)
    
    # append the predicted class
    predicted.cell.dat <- data.table(unlabelled.dat, pr)
    predicted.label <- paste0("Predicted_", label.col)
    
    names(predicted.cell.dat)[c(ncol(predicted.cell.dat))] <- predicted.label
    
    if (suppressWarnings(all(!is.na(as.numeric(as.character(predicted.cell.dat[[predicted.label]])))))) {
      predicted.cell.dat[[predicted.label]] <- as.numeric(as.character(predicted.cell.dat[[predicted.label]]))
    } 
    
    # don't just use as.numeric as it will mess up the order.
    #predicted.cell.dat[, ('Prediction') := as.numeric(levels(predicted.label))[predicted.label]]
    
    return(predicted.cell.dat)
  }
}
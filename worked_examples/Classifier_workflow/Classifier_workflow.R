##########################################################################################################
#### DRAFT Cytometry Analysis Pipeline for large and compleX data (CAPX) v3.0 - using the Spectre R package
#### Part 2 - Classification, save files
##########################################################################################################

# Givanna Putri, Thomas Myles Ashhurst
# 2020-05-18
# Spectre R package: https://sydneycytometry.org.au/spectre

##########################################################################################################
#### 1. Install packages, load packages, and set working directory
##########################################################################################################

library(Spectre)
package.check()
package.load()

library(caret)

dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

##########################################################################################################
#### Load data and run some clustering/UMAP etc
##########################################################################################################

cell.dat <- as.data.table(demo.start)
cell.dat <- do.subsample(cell.dat, 10000)

as.matrix(names(cell.dat))
CellularCols <- names(cell.dat)[c(2:10)]
CellularCols

as.matrix(names(cell.dat))
ClusterCols <- names(cell.dat)[c(2:10)]
ClusterCols

cell.dat <- run.flowsom(cell.dat, use.cols = ClusterCols, meta.k = 10)


## Check the data to see what it looks like.
cell.dat <- run.umap(cell.dat, use.cols = ClusterCols)

make.colour.plot(cell.dat, 
                 x.axis = "UMAP_X", 
                 y.axis = "UMAP_Y", 
                 col.axis = "FlowSOM_metacluster", 
                 col.type = 'factor',
                 add.label = TRUE,
                 save.to.disk = TRUE)


make.colour.plot(cell.dat, 
                 x.axis = "UMAP_X", 
                 y.axis = "UMAP_Y", 
                 col.axis = "FileName", 
                 col.type = 'factor',
                 add.label = TRUE,
                 save.to.disk = FALSE)

make.multi.plot(cell.dat, 
                x.axis = "UMAP_X", 
                y.axis = "UMAP_Y", 
                plot.by = CellularCols, 
                figure.title = "CellularCols",
                add.density = TRUE,
                save.each.plot = FALSE)


###############################################################################################
#### Train Classifier
###############################################################################################

# train classifier using various number of neighbours.
knn.stats <- Spectre::run.train.knn.classifier(dat = cell.dat,
                                               use.cols = ClusterCols,
                                               label.col = "FlowSOM_metacluster")

# see the training results and choose the suitable number of k.
knn.stats

###############################################################################################
#### Run classifier
###############################################################################################
## Predict the unlabelled data.
unpredicted.data <- as.data.table(demo.start)
unpredicted.data <- do.subsample(unpredicted.data, 10000)
predicted.data <- Spectre::run.knn.classifier(train.dat = cell.dat,
                                              unlabelled.dat = unpredicted.data,
                                              use.cols = ClusterCols,
                                              label.col = "FlowSOM_metacluster",
                                              num.neighbours = 3)

###############################################################################################
#### Assess classification
###############################################################################################

predicted.data <- run.umap(predicted.data, use.cols = ClusterCols)
make.colour.plot(predicted.data, 
                 x.axis = "UMAP_X", 
                 y.axis = "UMAP_Y",
                 col.type = "factor",
                 col.axis = "Prediction", 
                 add.label = TRUE,
                 save.to.disk = TRUE)


# Save data
Spectre::write.files(dat = predicted.data,
                     file.prefix= "Predicted_data", # required
                     write.csv = TRUE,
                     write.fcs = FALSE)

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
cell.dat <- do.subsample(cell.dat, "random", targets = 10000)

as.matrix(names(cell.dat))
CellularCols <- names(cell.dat)[c(2,4:6,8:9,11:13,16:19,21:30,32)]
CellularCols

as.matrix(names(cell.dat))
ClusterCols <- names(cell.dat)[c(11,13,17,18,19,21,22,23,24,27,28,32)]
ClusterCols

cell.dat <- run.flowsom(cell.dat, use.cols = ClusterCols, meta.k = 10)
cell.dat <- run.umap(cell.dat, use.cols = ClusterCols)
cell.dat <- run.tsne(cell.dat, use.cols = ClusterCols)

make.factor.plot(cell.dat, x.axis = "UMAP_X", y.axis = "UMAP_Y", col.axis = "FlowSOM_metacluster", add.label = TRUE)
make.factor.plot(cell.dat, x.axis = "UMAP_X", y.axis = "UMAP_Y", col.axis = "Group")
make.multi.marker.plot(cell.dat, x.axis = "UMAP_X", y.axis = "UMAP_Y", plot.by = CellularCols, figure.title = "CellularCols")
make.multi.marker.plot(cell.dat, x.axis = "UMAP_X", y.axis = "UMAP_Y", plot.by = ClusterCols, figure.title = "ClusterCols")


##########################################################################################################
#### Extract 1/3 of the dataset to use as 'training' data for the classifier (For demonstration only!)
#### When using workflow, this is where you need to specify the training data.
##########################################################################################################

head(cell.dat)
nrow(cell.dat)

thrd <- nrow(cell.dat)/3

set.seed(42)
rows <- sample(nrow(cell.dat))
cell.dat <- cell.dat[rows, ]

train.dat <- cell.dat[c(1:(2*thrd)),]
train.dat

##########################################################################################################
#### Train Classifier
##########################################################################################################

# train classifier using various number of neighbours.
knn.stats <- Spectre::run.train.knn.classifier(dat = train.dat,
                                      use.cols = ClusterCols,
                                      label.col = "FlowSOM_metacluster")

# see the training results and choose the suitable number of k.
knn.stats

##########################################################################################################
#### Run classifier
##########################################################################################################
## Predict the unlabelled data.
predicted.data <- Spectre::run.knn.classifier(train.dat = train.dat,
                                     unlabelled.dat = cell.dat,
                                     use.cols = ClusterCols,
                                     label.col = "FlowSOM_metacluster")

##########################################################################################################
#### Assess classification
##########################################################################################################

make.factor.plot(predicted.data, x.axis = "UMAP_X", y.axis = "UMAP_Y", col.axis = "Prediction", add.label = TRUE)


# Save data
Spectre::write.files(dat = cell.dat,
                     file.prefix= "Clustered_data", # required
                     write.csv = TRUE,
                     write.fcs = FALSE)

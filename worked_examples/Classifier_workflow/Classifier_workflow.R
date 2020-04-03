##########################################################################################################
#### DRAFT Cytometry Analysis Pipeline for large and compleX data (CAPX) v3.0 - using the Spectre R package
#### Part 2 - Classification, save files
##########################################################################################################

# Givanna Putri, Thomas Myles Ashhurst
# 2020-04-03
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

cell.dat <- run.flowsom(cell.dat, clustering.cols = ClusterCols, meta.clust.name = "FlowSOM_metacluster", clust.name = "FlowSOM_cluster", meta.k = 10)
cell.dat <- run.umap(cell.dat, use.cols = ClusterCols)
cell.dat <- run.tsne(cell.dat, use.cols = ClusterCols)

make.factor.plot(cell.dat, x.axis = "UMAP_X", y.axis = "UMAP_Y", col.axis = "FlowSOM_metacluster", add.label = TRUE)
make.factor.plot(cell.dat, x.axis = "UMAP_X", y.axis = "UMAP_Y", col.axis = "Group")
make.multi.marker.plot(cell.dat, x.axis = "UMAP_X", y.axis = "UMAP_Y", plot.by = CellularCols, figure.title = "CellularCols")
make.multi.marker.plot(cell.dat, x.axis = "UMAP_X", y.axis = "UMAP_Y", plot.by = ClusterCols, figure.title = "ClusterCols")


##########################################################################################################
#### Extract 1/3 of the dataset to use as 'training' data for the classifier
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

## Split the cellular markers columns and the identity of each cell into separate dataframe and vector
train.data.features <- train.dat[, ..ClusterCols]

# the label must be character label.
train.data.label <- as.character(train.dat[["FlowSOM_metacluster"]])

# train classifier using various number of neighbours
knn.stats <- train.knn.classifier(train.data = train.data.features, label = train.data.label)

knn.stats

##########################################################################################################
#### Run classifier
##########################################################################################################

## Now we use the classifier to predict the cell type of new data
unlabelled.data.features <- cell.dat[, ..ClusterCols]

predicted.label <- run.knn.classifier(train.data = train.data.features, # issue here -- the prediction labes come out in a character order, and are somehow paired on the basis of that order -- so cells in cluster '2' and assigned to predicted label '10' (which comes after 1, alphabetically)
                                      train.label = train.data.label,
                                      unlabelled.data = unlabelled.data.features,
                                      num.neighbours = 1)
new.dat <- cell.dat

# don't just use as.numeric as it will mess up the order.
new.dat[, ('Prediction') := as.numeric(levels(predicted.label))[predicted.label]]


##########################################################################################################
#### Assess classification
##########################################################################################################

make.factor.plot(cell.dat, x.axis = "UMAP_X", y.axis = "UMAP_Y", col.axis = "Prediction", add.label = TRUE)


# Save data
Spectre::write.files(dat = cell.dat,
                     file.prefix= "Clustered_data", # required
                     write.csv = TRUE,
                     write.fcs = FALSE)

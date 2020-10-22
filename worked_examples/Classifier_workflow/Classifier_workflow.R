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

dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# Load training data ====

cell.dat <- as.data.table(demo.clustered)
cell.dat <- do.subsample(cell.dat, 10000)

as.matrix(names(cell.dat))
classify_cols <- names(cell.dat)[c(11:19)]
classify_cols
label_col <- 'Population'

# Train classifier using various number of neighbours ====
knn.stats <- Spectre::run.train.knn.classifier(dat = cell.dat,
                                               use.cols = classify_cols,
                                               label.col = label_col)

# see the training results and choose the suitable number of k.
knn.stats

# Load un-labelled data
unpredicted.data <- as.data.table(demo.asinh)
unpredicted.data <- do.subsample(unpredicted.data, 1000)
predicted.data <- Spectre::run.knn.classifier(train.dat = cell.dat,
                                              unlabelled.dat = unpredicted.data,
                                              use.cols = classify_cols,
                                              label.col = label_col,
                                              num.neighbours = 9)
# Save data
Spectre::write.files(dat = predicted.data,
                     file.prefix= "Predicted_data", # required
                     write.csv = TRUE,
                     write.fcs = FALSE)

##########################################################################################################
#### DRAFT Cytometry Analysis Pipeline for large and compleX data (CAPX) v3.0 - using the Spectre R package
#### Part 2 - Classification, save files
##########################################################################################################

# Givanna Putri
# 2019-12-05
# Spectre R package: https://sydneycytometry.org.au/spectre

##########################################################################################################
#### 1. Install packages, load packages, and set working directory
##########################################################################################################

### 1.1. Load 'Spectre' package (using devtools)
if(!require('devtools')) {install.packages('devtools')}
library('devtools')

if(!require('Spectre')) {install_github("sydneycytometry/spectre", ref = "adding_classifier")}
library("Spectre")

### 1.2. Install packages

if(!require('class')) {install.packages('class')}
if(!require('caret')) {install.packages('caret')}
if(!require('rstudioapi')) {install.packages('rstudioapi')}
if(!require('data.table')) {install.packages('data.table')}

### 1.3. Load packages from library
library("class")
library("caret")
library("rstudioapi")
library("data.table")

### 1.4. Set working directory

## Set working directory
dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

## Can set manually using these lines, if desired
#PrimaryDirectory <- "/Users/thomasashhurst/Documents/Github/Public Github/Spectre/Other/Demo_dataset/"
#setwd(PrimaryDirectory)

## Create output directory
dir.create("Output_Classifier", showWarnings = FALSE)
setwd("Output_Classifier")
OutputDirectory <- getwd()
setwd(PrimaryDirectory)

##########################################################################################################
#### 2. Read and prepare data
##########################################################################################################
### Read SAMPLES (data) into workspace and review

## Start with reading in the labelled data i.e. the data to be used to train the classifier
TrainingDataDirectory <- paste(PrimaryDirectory, "training_data", sep="/")
setwd(TrainingDataDirectory)

## List of CSV files in PrimaryDirectory ## ADD A PRE-PROCESSING SCRIPT BEFORE THIS ONE -- FILE MERGE etc             ## HERE WE WANT ONE FILE PER SAMPLE
list.files(TrainingDataDirectory, ".csv")

## Import samples (read files into R from disk)
Spectre::read.files(file.loc = TrainingDataDirectory,
                    file.type = ".csv",
                    do.embed.file.names = TRUE)

ncol.check    # Review number of columns (features, markers) in each sample
nrow.check    # Review number of rows (cells) in each sample
name.table    # Review column names and their subsequent values

## Check data
head(data.list)
head(data.list[[1]])

## Save starting data
data.start <- data.list

## Merge files and review
Spectre::file.merge(x = data.list)

## Save as training data
train.data <- cell.dat

## Start with reading in the unlabelled data
UnlabelledDataDirectory <- paste(PrimaryDirectory, "unlabelled_data", sep="/")
setwd(UnlabelledDataDirectory)

## List of CSV files in PrimaryDirectory ## ADD A PRE-PROCESSING SCRIPT BEFORE THIS ONE -- FILE MERGE etc             ## HERE WE WANT ONE FILE PER SAMPLE
list.files(UnlabelledDataDirectory, ".csv")

## Import samples (read files into R from disk)
Spectre::read.files(file.loc = UnlabelledDataDirectory,
                    file.type = ".csv",
                    do.embed.file.names = TRUE)

ncol.check    # Review number of columns (features, markers) in each sample
nrow.check    # Review number of rows (cells) in each sample
name.table    # Review column names and their subsequent values

## Check data
head(data.list)
head(data.list[[1]])

## Save starting data
data.start <- data.list

### Merge files

## Merge files and review
Spectre::file.merge(x = data.list)

## Save as training data
unlabelled.data <- cell.dat


## Define key columns that might be used or dividing data (samples, groups, batches, etc)

file.col <- "Filename"
sample.col <- "Sample"
group.col <- "Group"
batch.col <- "Batch"

## Create a list of column names
ColumnNames <- as.matrix(unname(colnames(train.data))) # assign reporter and marker names (column names) to 'ColumnNames'
ColumnNames

### Define cellular and clustering columns

## Define columns that are 'valid' cellular markers (i.e. not live/dead, blank channels etc)
## To be used by the classifier to identify cell
ValidCellularColsNos <- c(5,6,8,9,11:13,16:19,21:30,32)
ValidCellularCols <- ColumnNames[ValidCellularColsNos]

ValidCellularCols  # check that the column names that appear are the ones you want to analyse

## Define the column that identify the population a cell belong to (i.e. Neutrophils, B cells)
## Note this column must NOT be numeric!
CellNameColNos <- c(39)
CellNameCol <- ColumnNames[CellNameColNos]

CellNameCol

## Split the cellular markers columns and the identity of each cell into separate dataframe and vector 
train.data.features <- train.data[, ValidCellularCols]
# the label must be character label.
train.data.label <- as.character(train.data[, CellNameCol])

# train classifier using various number of neighbours
knn.stats <- train.classifier(train.data = train.data.features, label = train.data.label)
knn.stats

## Now we use the classifier to predict the cell type of new data
unlabelled.data.features <- unlabelled.data[, ValidCellularCols]

predicted.label <- run.classifier(train.data = train.data.features,
               train.label = train.data.label,
               unlabelled.data = unlabelled.data.features,
               num.neighbours = 1)

unlabelled.data['Population'] <- predicted.label

# Save data
setwd(OutputDirectory)

Spectre::write.files(x = unlabelled.data,
                     file.prefix= "Clustered_data", # required
                     write.csv = TRUE,
                     write.fcs = FALSE)

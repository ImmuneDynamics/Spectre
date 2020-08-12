##########################################################################################################
#### Spectre -- Time Series Workflow
#### Part 1/1 - Clustering
##########################################################################################################

# Spectre R package: https://sydneycytometry.org.au/spectre
# Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### 1. Install packages, load packages, and set working directory
##########################################################################################################

library(Spectre)
Spectre::package.check()
Spectre::package.load() # If you do not have anaconda installed in your system, reticulate will prompt you to install one when you run this function.

## Set working directory
dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

## Can set manually using these lines, if desired
#PrimaryDirectory <- "/Users/Tom/Desktop/TAXXX"
#setwd(PrimaryDirectory)

## Create output directory
dir.create("Output_Spectre_cc", showWarnings = FALSE)
setwd("Output_Spectre_cc")
OutputDirectory <- getwd()
setwd(PrimaryDirectory)

##########################################################################################################
#### 2. Read and prepare data
##########################################################################################################

### Read SAMPLES (data) into workspace and review

## List of CSV files in input directory
InputDirectory <- paste0(PrimaryDirectory, "/data")
list.files(InputDirectory, ".csv")

## Import samples (read files into R from disk)
data.list <- Spectre::read.files(file.loc = InputDirectory,
                                 file.type = ".csv",
                                 do.embed.file.names = TRUE)

head(data.list)
head(data.list[[1]])

## Read in metadata
MetaDirectory <- paste0(PrimaryDirectory, '/metadata')
setwd(MetaDirectory)
meta.dat <- read.csv(file = "sample.details.csv")
meta.dat
setwd(PrimaryDirectory)

cell.dat <- rbindlist(data.list, fill = TRUE)
head(cell.dat)

## Embed sample metadata
to.add <- meta.dat[,c(1:4)]
to.add
cell.dat <- do.add.cols(cell.dat, "FileName", to.add, "Filename", rmv.ext = TRUE)

head(cell.dat)

# There should be 1 file per time point. IF you have 2 time points, you should have 2 items printed out!
# In this example, the time point is indicated by the Sample column.
as.matrix(unique(cell.dat[["Sample"]]))

##########################################################################################################
#### 3. Define data and sample variables for analysis
##########################################################################################################

### Define key columns

as.matrix(names(cell.dat))

## Define key columns that might be used or dividing data (samples, groups, batches, etc)
exp.name <- "TimeSeriesDemo"

timepoint.col <- "Sample"


## Create a list of column names
ColumnNames <- as.matrix(unname(colnames(cell.dat))) # assign reporter and marker names (column names) to 'ColumnNames'
ColumnNames

### Define columns for clustering

## Define columns that are 'valid' cellular markers (i.e. not live/dead, blank channels etc)
ClusteringColNos <- c(1:14)

ClusteringCols <- ColumnNames[ClusteringColNos] # e.g. [c(11, 23, 10)] to include the markers corresponding to the column numbers 11, 23, 10

ClusteringCols  # check that the column names that appear are the ones you want to analyse
ColumnNames[-ClusteringColNos] # Check which columns are being EXCLUDED!

##########################################################################################################
#### 4. Perform clustering
##########################################################################################################
## Prepare the environment where ChronoClust will execute.
# EXPERT USERS ONLY: if your anaconda is installed in custom location,
# pass the location to the function using parameter environment_path, like:
# environment_path = "/opt/conda/bin"
# NOTE: when running this in Docker, make sure you pass environment_path = "/opt/conda/bin"
# in addition to the 3 parameters.
Spectre::run.prepare.chronoclust(environment_name = "r-chronoclust",
                                 create_environment = FALSE,
                                 install_dependencies = FALSE)

# Set config for Chronoclust.
# We leave out k, lambda, pi, omicron, upsilon to default value.
config <- list(beta= 0.0,
               delta= 1.0,
               epsilon= 0.3,
               lambda= 0,
               k= 1,
               mu= 0.001,
               pi= 0,
               omicron= 0.0,
               upsilon= 1)

## Subsample for testing.
cell.dat.full <- cell.dat
cell.dat <- Spectre::do.subsample(dat = cell.dat,
                                      targets = 10000)

# if the following fail with the following error message:
# Error in py_module_import(module, convert = convert) : ModuleNotFoundError: No module named
# please manually restart your R-session and try again.
cell.dat <- Spectre::run.chronoclust(dat=cell.dat,
                                     timepoint.col=timepoint.col,
                                     use.cols=ClusteringCols,
                                     config=config)

# Check data
head(cell.dat)

# Write results out
setwd(OutputDirectory)
Spectre::write.files(cell.dat, exp.name)

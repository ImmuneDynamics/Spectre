##########################################################################################################
#### Spectre -- Time Series Workflow
#### Part 1/1 - Clustering
##########################################################################################################

# Spectre R package: https://sydneycytometry.org.au/spectre
# Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### 1. Install packages, load packages, and set working directory
##########################################################################################################

### 1.1. Install 'Spectre' package (using devtools) and the dependencies that Spectre requires

# For instructions on installing Spectre, please visit https://wiki.centenary.org.au/display/SPECTRE

### 1.2. Load packages

library(Spectre)
Spectre::package.check() # --> change so that message at the end is "All required packages have been successfully installed"
Spectre::package.load() # --> change so that message at the end is "All required packages have been successfully loaded"

session_info()

### 1.3. Set number of threads for data.table functions

getDTthreads()

### 1.4. Set working directory

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

## List of CSV files in PrimaryDirectory # HERE WE WANT ONE FILE PER SAMPLE
list.files(PrimaryDirectory, ".csv")

## Import samples (read files into R from disk)
data.list <- Spectre::read.files(file.loc = PrimaryDirectory,
                                 file.type = ".csv",
                                 do.embed.file.names = TRUE)

## Some checks
ncol.check    # Review number of columns (features, markers) in each sample
nrow.check    # Review number of rows (cells) in each sample
name.table    # Review column names and their subsequent values

head(data.list)
head(data.list[[1]])

### Merge files

## Merge files and review
cell.dat <- Spectre::do.merge.files(dat = data.list)

str(cell.dat)
head(cell.dat)
dim(cell.dat)

# There should be 1 file per time point. IF you have 2 time points, you should have 2 items printed out!
as.matrix(unique(cell.dat[["FileName"]]))

## Are there any NAs present in cell.dat? Yes if 'TRUE', no if 'FALSE'
any(is.na(cell.dat))

## Cleanup (not necessary, but recommended)
#rm(data.list, data.start, ncol.check, nrow.check, all.file.names, all.file.nums)

##########################################################################################################
#### 3. Define data and sample variables for analysis
##########################################################################################################

### Define key columns

as.matrix(names(cell.dat))

## Define key columns that might be used or dividing data (samples, groups, batches, etc)
exp.name <- "TimeSeriesDemo"

file.col <- "FileName"


## Create a list of column names
ColumnNames <- as.matrix(unname(colnames(cell.dat))) # assign reporter and marker names (column names) to 'ColumnNames'
ColumnNames

### Define columns for clustering

## Define columns that are 'valid' cellular markers (i.e. not live/dead, blank channels etc)
ColumnNames
# including SCA-1 as we're doing time series clustering
ClusteringColNos <- c(1:3)
ClusteringCols <- ColumnNames[ClusteringColNos] # e.g. [c(11, 23, 10)] to include the markers corresponding to the column numbers 11, 23, 10

ClusteringCols  # check that the column names that appear are the ones you want to analyse
ColumnNames[-ClusteringColNos] # Check which columns are being EXCLUDED!

### Checks

head(cell.dat)
ClusteringCols

## Export data containing only columns for clustering out to csv files
dir.create("Input_ChronoClust", showWarnings = FALSE)
setwd("Input_ChronoClust")
InputDirectoryChronoClust <- getwd()

# subset the data such that each time point is individual file.
# in this case, the column Time.point determine the time point
Spectre::do.split.data(cell.dat, 'FileName')


setwd(PrimaryDirectory)

##########################################################################################################
#### 4. Perform clustering
##########################################################################################################

## Get the input files for ChronoClust
chronoclust.input.files <- list.files(InputDirectoryChronoClust, ".csv")
chronoclust.input.files <- paste(InputDirectoryChronoClust, chronoclust.input.files, sep="/")

# Set config for Chronoclust.
# We leave out k, lambda, pi, omicron, upsilon to default value.
config <- list(beta= 0.2,
               delta= 0.05,
               epsilon= 0.03,
               mu= 0.01)


Spectre::run.prepare.chronoclust(environment_name = "chronoclust-R",
                            create_environment = TRUE,
                            install_dependencies = TRUE)


# if the following fail with the following error message:
# Error in py_module_import(module, convert = convert) : ModuleNotFoundError: No module named
# please manually restart your R-session and try again.
Spectre::run.chronoclust(data.files=chronoclust.input.files,
                output.dir=OutputDirectory,
                config=config)


# Assuming you have setup a conda environment, specify what's the environment name is
#conda.env.name <- "base"
#conda.env.location <- "/opt/conda/bin"

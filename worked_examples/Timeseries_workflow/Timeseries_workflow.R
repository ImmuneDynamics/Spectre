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

primary.dir <- 'M:/givanna/spectre_paper/WNV CNS timecourse (channel value)/'
source(paste0(primary.dir, 'run.chronoclust.R'))
source(paste0(primary.dir, 'make.network.plot.R'))

## Create output directory
output.dir <- paste0(primary.dir, "/Output_Spectre_cc_5")
dir.create(output.dir, showWarnings = FALSE)

##########################################################################################################
#### 2. Read and prepare data
##########################################################################################################

### Read SAMPLES (data) into workspace and review

## List of CSV files in input directory
data.dir <- paste0(primary.dir, "/data")
list.files(data.dir, ".csv")

## Import samples (read files into R from disk)
data.list <- Spectre::read.files(file.loc = data.dir,
                                 file.type = ".csv",
                                 do.embed.file.names = TRUE)


head(data.list)
head(data.list[[1]])

## Read in metadata
meta.dir <- paste0(primary.dir, '/metadata')
setwd(meta.dir)
meta.dat <- read.csv(file = "sample.details.csv")
meta.dat
setwd(primary.dir)

cell.dat <- rbindlist(data.list, fill = TRUE)
head(cell.dat)

## Embed sample metadata
to.add <- meta.dat[,c(1:4)]
to.add
cell.dat <- do.add.cols(cell.dat, "FileName", to.add, "Filename", rmv.ext = TRUE)

head(cell.dat)

# There should be 1 file per time point. IF you have 2 time points, you should have 2 items printed out!
# In this example, the time point is indicated by the Group column.
as.matrix(unique(cell.dat[["Group"]]))

##########################################################################################################
#### 3. Define data and sample variables for analysis
##########################################################################################################

### Define key columns

as.matrix(names(cell.dat))

## Define key columns that might be used or dividing data (samples, groups, batches, etc)
exp.name <- "WNV_brain"

timepoint.col <- "Group"
timepoints <- c("Mock", "WNV-01", "WNV-02", "WNV-03", "WNV-04", "WNV-05")

cluster.cols.nos <- c(1:8,10:20)
cluster.cols <- names(cell.dat)[cluster.cols.nos]

# not needed but just in case
cell.dat <- cell.dat[order(Group)]

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

# subsample the data if need be
# cell.dat[, .N, by=.(Group)]
cell.dat <- Spectre::do.subsample(cell.dat,
                                  rep(462803, 6),
                                  divide.by = 'Group')

# Set config for Chronoclust.
# We leave out k, lambda, pi, omicron, upsilon to default value.
config <- list(beta= 0.8,
               delta= 0.0,
               epsilon= 0.20,
               lambda= 0,
               k= 1,
               mu= 0.0005,
               pi= 0,
               omicron= 0.0,
               upsilon= 1)

setwd(output.dir)

# if the following fail with the following error message:
# Error in py_module_import(module, convert = convert) : ModuleNotFoundError: No module named
# please manually restart your R-session and try again.
cell.dat <- run.chronoclust(dat=cell.dat,
                            timepoint.col=timepoint.col,
                            timepoints=timepoints,
                            use.cols=cluster.cols,
                            config=config,
                            clean.up = FALSE)

# Check data
head(cell.dat)

unique(cell.dat$ChronoClust_cluster_lineage)

# Write results and parameters out
setwd(output.dir)
Spectre::write.files(cell.dat, exp.name)
sink(paste0("param_", exp.name, ".txt"))
print(config)
sink()

getwd()

list.files()

colnames(cell.dat) <- gsub(" ","_",colnames(cell.dat))
colnames(cell.dat) <- gsub("-","_",colnames(cell.dat))
timepoint.col <- "Group"
timepoints <- c("Mock", "WNV-01", "WNV-02", "WNV-03", "WNV-04", "WNV-05")
cluster.cols.nos <- c(1:8,10:20)

make.network.plot(dat = cell.dat,
                  timepoint.col = timepoint.col,
                  timepoints = timepoints,
                  cluster.col = 'ChronoClust_cluster_lineage',
                  marker.cols = cluster.cols.nos,
                  node.size = 6,
                  arrow.length = 2,
                  arrow.head.gap = 3)

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

### Set or create 'input' directory

setwd(primary.dir)
dir.create('../data', showWarnings = FALSE)
setwd("../data/")
InputDirectory <- getwd()
setwd(primary.dir)

### Set or create 'metadata' directory

setwd(primary.dir)
dir.create('../metadata', showWarnings = FALSE)
setwd("../metadata/")
MetaDirectory <- getwd()
setwd(primary.dir)

## Create output directory
setwd(primary.dir)
output.dir <- paste0(primary.dir, "/Output_Spectre_cc")
dir.create(output.dir, showWarnings = FALSE)

##########################################################################################################
#### 2. Read and prepare data
##########################################################################################################

### If you need the demo dataset, uncomment the following code (select all, CMD+SHIFT+C) and run to download
### Alternative: download from https://github.com/ImmuneDynamics/data/blob/main/timeSeries.zip?raw=TRUE

    # setwd(primary.dir)
    # setwd("../")
    # getwd()
    # download.file(url = "https://github.com/ImmuneDynamics/data/blob/main/timeSeries.zip?raw=TRUE", destfile = 'timeSeries.zip', mode = 'wb')
    # unzip(zipfile = 'timeSeries.zip')
    # for(i in list.files('timeSeries/data', full.names = TRUE)){
    #   file.rename(from = i,  to = gsub('timeSeries/', '', i))
    # }
    # for(i in list.files('timeSeries/metadata', full.names = TRUE)){
    #   file.rename(from = i,  to = gsub('timeSeries/', '', i))
    # }
    # unlink(c('timeSeries/', 'timeSeries.zip', '__MACOSX'), recursive = TRUE)

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
# cell.dat <- Spectre::do.subsample(cell.dat,
#                                   rep(462803, 6),
#                                   divide.by = 'Group')

# compute the suitable epsilon. find the kink and convert to between 0 and 1.
# TODO: beta version. use and interpret with caution!
#dat.sub.mat <- as.matrix(sapply(dat.sub, as.numeric))
#dbscan::kNNdistplot(dat.sub.mat, k=1)

# Set config for Chronoclust.
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

#### Draw plot ####
cell.dat.plot <- cell.dat[cell.dat$ChronoClust_cluster_lineage != 'None',]

colnames(cell.dat.plot) <- gsub(" ","_",colnames(cell.dat.plot))
colnames(cell.dat.plot) <- gsub("-","_",colnames(cell.dat.plot))

setwd(output.dir)
dir.create("spectral")
setwd("spectral")
make.network.plot(dat = cell.dat.plot,
                  timepoint.col = timepoint.col,
                  timepoints = timepoints,
                  cluster.col = 'ChronoClust_cluster_lineage',
                  marker.cols = cluster.cols.nos,
                  node.size = 'auto',
                  arrow.length = 3,
                  arrow.head.gap = 2)

### Draw plot where the node size is the proportion of cells ###
setwd(primary.dir)
dir.create("inferno_static")
setwd("inferno_static")
make.network.plot(dat = cell.dat.plot,
                  timepoint.col = timepoint.col,
                  timepoints = timepoints,
                  cluster.col = 'ChronoClust_cluster_lineage',
                  marker.cols = cluster.cols.nos,
                  node.size = 6,
                  arrow.length = 3,
                  arrow.head.gap = 2,
                  standard.colours = 'inferno')

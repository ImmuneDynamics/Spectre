##################################################################################################
##### ChronoClust basic tutorial #####################################################################
##################################################################################################

# Givanna Putri
# 2019-12-06
# Spectre R package: https://sydneycytometry.org.au/spectre

##########################################################################################################
#### 1. Install packages, load packages
##########################################################################################################

library("Spectre")

Spectre::check.packages()
Spectre::load.packages()

##########################################################################################################
#### 2. Setup the working directory and the config files for Chronoclust
##########################################################################################################

## Set working directory
dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()

# Using test data stored in data_files directory
data.files <- c('sample_data_d0.csv',
                'sample_data_d1.csv')

# Set config for Chronoclust.
config <- list(beta= 0.2,
               delta= 0.05,
               epsilon= 0.03,
               lambda= 2,
               k= 4,
               mu= 0.01,
               pi= 3,
               omicron= 0.00000435,
               upsilon= 6.5)

# Set directory where Chronoclust will store results
outdir <- "output_chronoclust"

##########################################################################################################
#### 3. Run Chronoclust
##########################################################################################################

# assuming you have setup a conda environment, specify what's the environment name is
# conda.env.name <- "r-reticulate" # note, 'anaconda' (https://www.anaconda.com/distribution/#download-section) needs to be installed for the conda environment to work
conda.env.name <- "chronoclust-public"
# location where anaconda is installed in your computer. Read more below.
conda.env.location <- NULL
# OPTIONAL: conda is normally installed under /Users/<yourname>/anaconda3, but it can also be installed
# elsewhere such as /opt/anaconda3/bin/conda.
# Hence if you encounter error such as "Have you installed conda?", then it's likely that your conda is not
# installed at the default location. Thus you need to find where it's installed and specify it
# conda.env.location <- "/root/miniconda3/bin/conda"

Spectre::run.chronoclust(data.files=data.files,
                         output.dir=outdir,
                         conda.env.name=conda.env.name,
                         conda.env.location=conda.env.location,
                         config=config)



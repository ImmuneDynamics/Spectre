##########################################################################################################
#### DRAFT Cytometry Analysis Pipeline for large and compleX data (CAPX) v3.0 - using the Spectre R package
#### Part 3 - Clustering with Chronoclust, save files
##########################################################################################################

# Givanna Putri
# 2019-12-06
# Spectre R package: https://sydneycytometry.org.au/spectre

##########################################################################################################
#### 1. Install packages, load packages
##########################################################################################################

### 1.1. Load 'Spectre' package (using devtools)
if(!require('devtools')) {install.packages('devtools')}
library('devtools')

# Master version
#if(!require('Spectre')) {install_github("sydneycytometry/spectre")}

# Other branch version. "adding_chronoclust" is the branch name.
devtools::install_github(repo = "sydneycytometry/spectre", ref = 'adding_chronoclust')

library("Spectre")

### 1.2 Load 'Reticulate' package
if(!require('reticulate')) {install.packages('reticulate')}
library("reticulate")

##########################################################################################################
#### 2. Setup the working directory and the config files for Chronoclust
##########################################################################################################

setwd("/Users/givanna/Documents/phd/code/Spectre/worked_example/Chronoclust_workflow")
# specify the location of the xml file containing the dataset for each time point
input.xml.location <- "config_files/input.xml"

# specify the location of the xml file containing the parameter values for chronoclust to operate on
config.xml.location <- "config_files/config.xml"

# specify the output directory for chronoclust
out.dir <- "Output_chronoclust"

##########################################################################################################
#### 3. Run Chronoclust
##########################################################################################################
Spectre::run.chronoclust(config.xml.file.location = config.xml.location,
                         input.xml.file.location = input.xml.location,
                         output.dir.location = out.dir)




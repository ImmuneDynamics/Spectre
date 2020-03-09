## Install (if not already installed)
if(!require('devtools')) {install.packages('devtools')}

## Install Spectre
library('devtools')
install_github("sydneycytometry/spectre")

## Install BiocManager to download packages from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

## Download additional BioConductor packages
if(!require('flowCore')) {BiocManager::install('flowCore')}
if(!require('Biobase')) {BiocManager::install('Biobase')}
if(!require('flowViz')) {BiocManager::install('flowViz')}
if(!require('FlowSOM')) {BiocManager::install('FlowSOM')}

## Check if all required packages have been installed
Spectre::package.check()

## Load all required packages
Spectre::package.load()

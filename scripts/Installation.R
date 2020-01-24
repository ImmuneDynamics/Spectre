##########################################################################################################
#### Install Spectre and dependent packages
##########################################################################################################

### 1.1. Install Spectre and dependent packages
    
    ## Install devtolls
    if(!require('devtools')) {install.packages('devtools')}
    library('devtools')
    
    ## Install Spectre
    #install_github("sydneycytometry/spectre")
    install_github("sydneycytometry/spectre", ref = 'development') # option to install the development verison if required
    library("Spectre")
    
    ## Install BiocManager to download packages from Bioconductor
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    ## Download additional BioConductor packages
    if(!require('flowCore')) {BiocManager::install('flowCore')}
    if(!require('Biobase')) {BiocManager::install('Biobase')}
    if(!require('flowViz')) {BiocManager::install('flowViz')}
    if(!require('FlowSOM')) {BiocManager::install('FlowSOM')}

### 1.2. Load packages to ensure successful installation
    
    library(Spectre)
    Spectre::check.packages() # --> change so that message at the end is "All required packages have been successfully installed"
    Spectre::load.packages() # --> change so that message at the end is "All required packages have been successfully loaded"
    
    session_info()
    
    
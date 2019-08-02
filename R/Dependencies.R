

    ### Install and load packages

    ## From CRAN (required)
    if(!require('plyr')) {install.packages('plyr')}
    if(!require('data.table')) {install.packages('data.table')}
    #install.packages('data.table')

    if(!require('rstudioapi')) {install.packages('rstudioapi')}
    if(!require('devtools')){install.packages("devtools")}

    ## From CRAN (required for clustering and tSNE)
    if(!require('FlowSOM')) {source("https://bioconductor.org/biocLite.R") # for running FlowSOM ### POSSIBLY INSTALL FROM GITHUB
      biocLite('FlowSOM')}
    if(!require('Rtsne')) {install.packages("Rtsne")} # for running tSNE

    ## From CRAN (only required for plotting -- if installation unsuccessful, set Run_tSNEplots and Run_ClusterPlots to 0)
    if (!require("ggplot2")){install.packages("ggplot2")} # for plotting tSNE graphs
    if (!require("colorRamps")){install.packages("colorRamps")} # for colour scheme management
    if (!require("ggthemes")){install.packages("ggthemes")} # for plot themes
    if (!require("scales")){install.packages("scales")} # for re-scaling if necessary

    ## From Bioconductor (required)
    if(!require('flowViz')) {source("https://bioconductor.org/biocLite.R")
      biocLite('flowViz')}
    if(!require('flowCore')) {source("https://bioconductor.org/biocLite.R")
      biocLite('flowCore')}
    if(!require('Biobase')) {source("https://bioconductor.org/biocLite.R")
      biocLite('Biobase')}

    ## From Github repositories
    #if(!require('Rphenograph')){devtools::install_github("JinmiaoChenLab/Rphenograph")} # not required

    ### Load packages
    library('plyr')
    library('data.table')
    library('rstudioapi')
    library('devtools')
    library('FlowSOM') ###
    library('Rtsne')
    library('ggplot2')
    library('colorRamps')
    library('ggthemes')
    library('scales')
    library('flowCore') ###
    library('Biobase')
    library('flowViz') ### matrixStats




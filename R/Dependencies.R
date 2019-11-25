### Install and load packages

    ## From CRAN (required)
    if(!require('plyr')) {install.packages('plyr', repos = "https://mirror.aarnet.edu.au/pub/CRAN/")}
    if(!require('data.table')) {install.packages('data.table', repos = "https://mirror.aarnet.edu.au/pub/CRAN/")}
    if(!require('tidyr')) {install.packages('tidyr', repos = "https://mirror.aarnet.edu.au/pub/CRAN/")}
    if(!require('rstudioapi')) {install.packages('rstudioapi')}
    #if(!require('rtsne')) {install.packages('rtsne')}

    install.packages('rtsne')

    ### Load packages
    library('plyr')
    library('data.table')
    library('tidyr')
    library('rstudioapi')
    library('rtsne')

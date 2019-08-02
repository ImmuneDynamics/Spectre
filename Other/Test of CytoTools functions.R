##### Performing UMAP on cytometry data using the Spectre package #####
# Thomas Myles Ashhurst
# 2019-08-02
# https://sydneycytometry.org.au/spectre


### 1. Install packages, load packages, and set workin directory

    ## 1.1. Install packages
    #if(!require('Spectre')) {install.packages('Spectre')}
    if(!require('plyr')) {install.packages('plyr')}
    if(!require('data.table')) {install.packages('data.table')}
    if(!require('rstudioapi')) {install.packages('rstudioapi')}
    if(!require('devtools')) {install.packages('devtools')}
    if(!require('umap')) {install.packages('umap')}

    ## 1.2. Load packages
    library('plyr')
    library('data.table')
    library('rstudioapi')
    library('devtools')
    library('umap')

    ## 1.3. Load 'Spectre' package
    if(!require('sydneycytometry/spectre')) {install_github("sydneycytometry/spectre")}
    library("Spectre")

    ## 1.4. Set working directory
    PrimaryDirectory <- "/Users/Tom/Downloads/CAPX-2.5/Demo dataset/"
    setwd(PrimaryDirectory)

    list.files(PrimaryDirectory, ".csv")

### 2. Read and prepare data

    ## Read in samples
        Spectre::read.files(file.loc = PrimaryDirectory, file.type = ".csv", do.embed.file.names = TRUE)

    ## Review sample properties
        ncol.check
        nrow.check

        name.table
        head(data.list)
        head(data.list[[1]])

    ## Add sample names/numbers

        #CytoTools::file.annotate()

    ## Add groups
        as.matrix(names(data.list))

            group.names = list()
            group.nums = list()

            group.names[[1]] <- "Mock"
            group.names[[2]] <- "WNV"

            group.nums[[1]] <- c(1:6)
            group.nums[[2]] <- c(7:12)

        Spectre::add.groups(x = data.list, grp.names = group.names, grp.nums = group.nums)


        head(data.list)
        head(data.list[[7]])

    ## Merge samples
    Spectre::merge.files(x = data.list)

    head(cell.dat)

    # Look for any NAs

    # Cleanup (not necessary, but recommended)
        #rm(data.list)
        #rm(group.names)
        #rm(group.nums)
        #rm(name.table)
        #rm(ncol.check)
        #rm(nrow.check)

### Perform FlowSOM clustering



## Downsample data in preparation for dimensionality reduction

    Spectre::subsample(x = cell.dat,
                       method = "per.sample",
                       samp.col = "FileName",
                       targets = c(rep(100, 12)),
                       seed = 42)

    Spectre::subsample(x = cell.dat,
                       method = "random",
                       samp.col = "FileName",
                       targets = c(100),
                       seed = 42)

    dim(subsample.res)

    cell.dat.sub <- subsample.res
    rm(subsample.res)




### 3. Perform UMAP

    as.matrix(names(cell.dat.sub))

    #x <- as.data.frame(cell.dat.sub)
    use.cols <- c(5,6,8,9)

    #Spectre::umap(x = cell.dat,
    #              use.cols = c(5,6,8,9,11,12,13,17:19,21:30,32),
    #              umap.seed = 42)



    #umap.to.plot <- cbind(cell.dat, umap.res)

    xdat <- cell.dat.sub
    xdat <- iris
    use.cols <- c(1:4)

    umap.seed <- 42

    xdat <- xdat[c(1:4)]

    custom.config <- umap.defaults
    custom.config$random_state <- umap.seed #umap.seed
    res <- umap(d = xdat[use.cols],
                config = custom.config)

    umap.res <- res$layout
    head(umap.res)

    umap.res <- as.data.frame(umap.res)
    head(umap.res)

    names(umap.res)[names(umap.res) == "V1"] <- paste0("UMAP", "_", umap.seed, "_", "X")
    names(umap.res)[names(umap.res) == "V2"] <- paste0("UMAP", "_", umap.seed, "_", "Y")

    head(umap.res)


    cell.dat.sub <- cbind(iris, umap.res)
    #plot.dat$umap

    head(cell.dat.sub)

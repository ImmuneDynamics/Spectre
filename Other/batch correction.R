##### Batch alignment using gaussNorm function #####

    # Thomas Myles Ashhurst
    # 2019-08-02
    # https://sydneycytometry.org.au

### 1. Install packages, load packages, and set working directory

    ## 1.1. Install packages, if they are not already installed
    if(!require('flowStats')) {install.packages('flowStats')}
    if(!require('flowCore')) {install.packages('flowCore')}
    if(!require('flowViz')) {install.packages('flowViz')}
    if(!require('Biobase')) {install.packages('Biobase')}

    ## 1.2. Load packages from library
    library(flowStats)
    library(flowCore)
    library(flowViz)
    library(Biobase)

    ## 1.3. Set working directory
    setwd("/Users/Tom/Downloads/CAPX-2.5/Demo dataset/Output_CSV-to-FCS_2019-08-02_16-20-05/")
    PrimaryDirectory <- getwd()

### 2. Read FCS files into flowSet

    ## 2.1. Generate a list of filenames for FCS files
    file.names <- list.files(path = PrimaryDirectory, pattern = ".fcs")

    ## 2.2. Read FCS files into a 'flowSet'
    dat.fs <- read.flowSet(files = file.names, path = PrimaryDirectory)

    ## 2.3. Review what is in the flowSet
    dat.fs
    ## 2.4. Transform axis
    #dat.fs <- transform(dat.fs,
    #                    "CD4"=asinh(CD4),
    #                    "CD3"=asinh(CD3),
    #                    "CD8"=asinh(CD8))

### 3. Perform gaussNorm batch correction

    ## 3.1. Review the column names in the flowSet
    as.matrix(dat.fs@colnames)

    ## 3.2. Choose which columns will be batch-corrected
    to.correct.num <- c(5:32)
    to.correct.names <- dat.fs@colnames[to.correct.num]

    ## 3.3. Perform gaussNorm correction on selected columns
    res <- gaussNorm(dat.fs, to.correct.names, max.lms = 1, peak.density.thr = 0.5)$flowset


### 4. Review results

    ## 4.1. Review column names
    as.matrix(dat.fs@colnames)

    ## 4.2. Choose the column to review 'original' and 'batch corrected' versions
    d1 <- densityplot(~BUV615.Siglec.F, dat.fs, main="original", filter=curv1Filter("BUV615.Siglec.F"))
    d2 <- densityplot(~BUV615.Siglec.F, res, main="batch-corrected", filter=curv1Filter("BUV615.Siglec.F"))

    ## 4.3. Plot the 'original' and 'batch-correced' plots for the selected column
    plot(d1, split=c(1,1,2,1))
    plot(d2, split=c(2,1,2,1), newpage=FALSE)


    ## Add a loop to create and save each image


### 5.Write batch-corrected FCS files

    ## 5.1. Create an output folder
    setwd(PrimaryDirectory)
    dir.create("Output", showWarnings = FALSE)
    setwd("Output")
    OutputDirectory <- getwd()

    ## 5.2. Now write out into  files
    write.flowSet(x = res, outdir = OutputDirectory)


    dat.fs[[2]]

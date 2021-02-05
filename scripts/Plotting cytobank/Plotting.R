##########################################################################################################
#### Spectre - custom plotting
##########################################################################################################

    # Spectre R package: https://sydneycytometry.org.au/spectre
    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### 1. Load packages, and set working directory
##########################################################################################################

    ### Load libraries

        library(Spectre)
        Spectre::package.check()    # Check that all required packages are installed
        Spectre::package.load()     # Load required packages

    ### Set PrimaryDirectory
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()

    ### Set 'input' directory
        setwd(PrimaryDirectory)
        InputDirectory <- getwd()
        setwd(PrimaryDirectory)

    ### Create output directory
        dir.create("Output_Spectre", showWarnings = FALSE)
        setwd("Output_Spectre")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### 2. Import and prep data
##########################################################################################################

    ### Import data

        setwd(InputDirectory)
        list.files(InputDirectory, ".fcs")

        data.list <- Spectre::read.files(file.loc = InputDirectory,
                                         file.type = ".fcs",
                                         do.embed.file.names = TRUE)

    ### Check the data

        check <- do.list.summary(data.list)

        check$name.table # Review column names and their subsequent values
        check$ncol.check # Review number of columns (features, markers) in each sample
        check$nrow.check # Review number of rows (cells) in each sample

        data.list[[1]]

    ### Merge data

        cell.dat <- Spectre::do.merge.files(dat = data.list)
        cell.dat
        
##########################################################################################################
#### 3. Data transformation
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 1 - transformed plots")
    setwd("Output 1 - transformed plots")
        
    ### Arcsinh transformation

        as.matrix(names(cell.dat))

        to.asinh <- names(cell.dat)[c(1:9)]
        to.asinh

        cofactor <- 500
        plot.against <- "Ly6C_asinh"

        cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
        transformed.cols <- paste0(to.asinh, "_asinh")

        for(i in transformed.cols){
          make.colour.plot(do.subsample(cell.dat, 1000), i, plot.against)
        }

##########################################################################################################
#### Re-scaling between 0 and 1
##########################################################################################################
        
    ### Re-scaling preferences
        
        cell.dat
        
        as.matrix(names(cell.dat))
        to.rescale <- names(cell.dat)[c(14:22)] # Choose the columns to rescale
        to.rescale
        
    ### Perform and check re-scaling
        
        cell.dat <- do.rescale(cell.dat, to.rescale)
        
    ### Perform and check re-scaling
        
        as.matrix(names(cell.dat))
        cell.dat

##########################################################################################################
#### Add group data
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 2 - plots")
    setwd("Output 2 - plots")
        
    ### Create multi plots
        
        as.matrix(names(cell.dat))
        
        plot.cols <- names(cell.dat)[c(14:22)] 
        plot.cols

        make.multi.plot(dat = cell.dat, x.axis = "FItSNE_X", y.axis = "FItSNE_Y", plot.by = plot.cols) # Change X and Y axis to desired tSNE columns
        
        for(i in plot.cols){
            make.multi.plot(dat = cell.dat, x.axis = "FItSNE_X", y.axis = "FItSNE_Y", plot.by = i, divide.by = "FileName") # Change X and Y axis to desired tSNE columns
        }

        

##########################################################################################################
#### Spectre: General Discovery Workflow - (1/4) - Data preparation and arcsinh transformation
##########################################################################################################

    # Spectre R package: https://github.com/ImmuneDynamics/Spectre
    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### Analysis session setup
##########################################################################################################

    ### Load packages

        library(Spectre)
        Spectre::package.check()    # Check that all required packages are installed
        Spectre::package.load()     # Load required packages

    ### Set DT threads

        getDTthreads()

    ### Set primary directory

        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

    ### Set input directory

        setwd(PrimaryDirectory)
        setwd("data/")
        InputDirectory <- getwd()
        InputDirectory
        setwd(PrimaryDirectory)

    ### Set metadata directory

        setwd(PrimaryDirectory)
        setwd("metadata/")
        MetaDirectory <- getwd()
        MetaDirectory
        setwd(PrimaryDirectory)

    ### Set output directory

        setwd(PrimaryDirectory)
        dir.create("Output 1 - data prep", showWarnings = FALSE)
        setwd("Output 1 - data prep")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### Import data
##########################################################################################################
    
    setwd(InputDirectory)
        
    ### Import data
        
        list.files(InputDirectory, ".csv")
        
        data.list <- Spectre::read.files(file.loc = InputDirectory,
                                         file.type = ".csv",
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
#### Acrsinh transformations
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 1 - transformation plots")
    setwd("Output 1 - transformation plots")

    ### Co-factor targets
    
        cf.low <- 100
        cf.mid <- 500
        cf.high <- 1000
        
        cf.low
        cf.mid
        cf.high
    
    ### Transformation settings

        as.matrix(names(cell.dat))

        asinh.low <- names(cell.dat)[c(1:3)]
        asinh.mid <- names(cell.dat)[c(4:6)]
        asinh.high <- names(cell.dat)[c(7:9)]
    
        asinh.low
        asinh.mid
        asinh.high
            
    ### Test transformation settings on subsampled data
        
        sub <- do.subsample(cell.dat, 10000)

        sub <- do.asinh(sub, use.cols = asinh.low, cofactor = cf.low)
        sub <- do.asinh(sub, use.cols = asinh.mid, cofactor = cf.mid)
        sub <- do.asinh(sub, use.cols = asinh.high, cofactor = cf.high)
        
        as.matrix(names(sub))

        transf.cols <- names(sub)[grepl('_asinh', names(sub))]
        transf.cols

    ### Make plots of transformed columns from the subsampled data
        
        plot.against <- 'Ly6C_asinh'
        which(names(sub) == plot.against)
        
        for(i in transf.cols){
          make.colour.plot(sub, i, plot.against)
        }
        
    ### Apply transformation to full dataset

        cell.dat <- do.asinh(cell.dat, use.cols = asinh.low, cofactor = cf.low)
        cell.dat <- do.asinh(cell.dat, use.cols = asinh.mid, cofactor = cf.mid)
        cell.dat <- do.asinh(cell.dat, use.cols = asinh.high, cofactor = cf.high)
        
        as.matrix(names(cell.dat))

##########################################################################################################
#### Add metadata and set preferences
##########################################################################################################
        
    setwd(MetaDirectory)
        
    ### Read in sample metadata
        
        meta.dat <- fread("sample.details.csv")
        meta.dat
        
        sample.info <- meta.dat[,c(1:4)]
        sample.info

    ### Add sample metadata to primary data.table
        
        cell.dat <- do.add.cols(cell.dat, "FileName", sample.info, "Filename", rmv.ext = TRUE)
        cell.dat

##########################################################################################################
#### Write data to disk
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 2 - transformed data")
    setwd("Output 2 - transformed data")
        
    ### Write cellular data and analysis  preferences to disk

        fwrite(cell.dat, "cell.dat.csv") # data

    ### Save session info to disk

        setwd(OutputDirectory)
        dir.create("Output - info", showWarnings = FALSE)
        setwd("Output - info")

        sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
        session_info()
        sink()

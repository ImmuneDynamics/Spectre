##########################################################################################################
#### Spectre discovery workflow - A - data preparation
##########################################################################################################

    # Spectre R package: https://sydneycytometry.org.au/spectre
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
        dir.create("Output A - data prep", showWarnings = FALSE)
        setwd("Output A - data prep")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### Import data
##########################################################################################################

    setwd(InputDirectory)

    ### Import data files

        data.list <- list()

        files <- list.files(getwd(), ".csv")
        files

        for(i in files){
          data.list[[i]] <- fread(i)
          nms <- names(data.list[[i]])

          data.list[[i]]$FileName <- i
          data.list[[i]] <- data.list[[i]][,c("FileName", nms), with = FALSE]
        }

        data.list[[1]]

    ### Check data and merge into single data.table

        check <- do.list.summary(data.list)

        check$name.table
        check$ncol.check
        check$nrow.check

        cell.dat <- rbindlist(data.list, fill = TRUE)
        cell.dat

##########################################################################################################
#### Transformations
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output - transformation plots")
    setwd("Output - transformation plots")

    ### Test transformations

        as.matrix(names(cell.dat))

        to.norm <- names(cell.dat)[c()]
        to.asinh <- names(cell.dat)[c(2:10)]

        sub <- do.subsample(cell.dat, 10000)

        sub <- do.normalise(sub, use.cols = to.norm, new.min = 0, new.max = 5)
        sub

        sub <- do.asinh(sub, use.cols = to.asinh, cofactor = 1000)
        sub

        as.matrix(names(sub))

        transf.cols <- names(sub)[c(11:19)]
        transf.cols

    ### Make plots

        for(i in transf.cols){
          make.colour.plot(sub, i, "Ly6C_asinh")
        }

    ### Apply to sample

        cell.dat <- do.normalise(cell.dat, use.cols = to.norm, new.min = 0, new.max = 5)
        cell.dat <- do.asinh(cell.dat, use.cols = to.asinh, cofactor = 1000)

        cell.dat
        as.matrix(names(cell.dat))

##########################################################################################################
#### Add metadata and set preferences
##########################################################################################################

    setwd(OutputDirectory)
    prefs <- list()

    ### Add critical metadata

        setwd(MetaDirectory)

        sample.dat <- fread("sample.details.csv")
        prefs$sample.dat <- sample.dat

        as.matrix(names(sample.dat))

        to.add <- sample.dat[,c(1:4)]
        to.add

        cell.dat <- do.add.cols(cell.dat, "FileName", to.add, "Filename")
        cell.dat

    ### Sample preferences

        as.matrix(names(cell.dat))

        prefs$sample.col <- "Sample"
        prefs$group.col <- "Group"
        prefs$batch.col <- "Batch"

    ### Downsample preferences
        as.matrix(unique(cell.dat$Group))

        prefs$sub.by <- "Group"
        prefs$sub.targets <- c(2000, 20000)

    ### Clustering preferences

        ## Cellular cols
        as.matrix(names(cell.dat))

        prefs$cellular.cols <- names(cell.dat)[c(11:19)]
        prefs$cellular.cols

        ## Columns for clustering
        as.matrix(names(cell.dat))

        prefs$clustering.cols <- names(cell.dat)[c(11:19)]
        prefs$clustering.cols

        ## Cluster numbers etc
        prefs$metak <- 'auto'

##########################################################################################################
#### Write data to disk
##########################################################################################################

    ### Write cellular data and analysis  preferences to disk

        setwd(OutputDirectory)

        fwrite(cell.dat, "Cellular.data.csv") # data
        saveRDS(prefs, "Analysis preferences.rds") # analysis preferences

    ### Save session info to disk

        setwd(OutputDirectory)
        dir.create("Output - info", showWarnings = FALSE)
        setwd("Output - info")

        sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
        session_info()
        sink()

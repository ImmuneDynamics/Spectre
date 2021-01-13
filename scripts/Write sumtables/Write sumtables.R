##########################################################################################################
#### Spectre - write sumtables
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
        InputDirectory <- getwd()
        InputDirectory
        setwd(PrimaryDirectory)

    ### Set output directory

        setwd(PrimaryDirectory)
        dir.create("Write sumtables", showWarnings = FALSE)
        setwd("Write sumtables")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### Import data
##########################################################################################################

    ### Data input
        # dat <- Spectre::demo.clustered
        dat <- fread()

    ### Specifications
        as.matrix(names(dat))


##########################################################################################################
#### Split and write
##########################################################################################################

    ###


    ###


    ###










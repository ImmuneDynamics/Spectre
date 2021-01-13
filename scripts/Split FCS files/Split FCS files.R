##########################################################################################################
#### Spectre - Split into FCS files
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
        dir.create("Split FCS files", showWarnings = FALSE)
        setwd("Split FCS files")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### Import data
########################################################################################################## 
    
    ### Read in data
        dat <- fread("AvC.csv")
        as.matrix(names(dat))
    
    ### Specify
        exp <- 'AvC'
        group.col <- 'DiseaseStatus'
        group.col2 <- 'DiseaseSeverity'
        sample.col <- 'Sample'

##########################################################################################################
#### Split and write
########################################################################################################## 
        
    ### Write FCS files
        setwd(OutputDirectory)
        
        for(i in unique(dat[[sample.col]])){
          # i <- unique(dat[[sample.col]])[[1]]
          temp <- dat[dat[[sample.col]] == i,]
          
          grp <- unique(temp[[group.col]])
          grp2 <- unique(temp[[group.col2]])
          samp <- unique(temp[[sample.col]])
          
          write.files(temp, file.prefix = paste0(exp, "_", grp, "_", grp2, "_", samp), write.csv = FALSE, write.fcs = TRUE)
        }

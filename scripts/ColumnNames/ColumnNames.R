##################################################################################
    
    ### Load libraries
    library(Spectre)
    Spectre::package.check()    # Check that all required packages are installed
    Spectre::package.load()     # Load required packages
    
    ### Set PrimaryDirectory
    dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
    getwd()
    PrimaryDirectory <- getwd()
    PrimaryDirectory
    
    ### Create output directory
    dir.create("Output - ColumnNames", showWarnings = FALSE)
    setwd("Output - ColumnNames")
    OutputDirectory <- getwd()
    OutputDirectory

##################################################################################

    setwd(PrimaryDirectory)
    
    files <- list.files(getwd(), '.csv')
    as.matrix(files)

    data.list <- list()
    
    for(i in files){
      data.list[[i]] <- fread(i, nrows = 1)
    }
    
    dat <- rbindlist(data.list, fill = TRUE)
    dat 
    
    as.matrix(names(dat))
    
##################################################################################    

    setwd(OutputDirectory)
    
    dt <- data.table('Column names' = names(dat), 'New names' = names(dat))     
    dt
    
    fwrite(dt, file = 'ColumnNames.csv')


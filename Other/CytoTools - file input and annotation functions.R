#################################################################################
###### CytoTools - file unput and anootation functions
#################################################################################

  # File input and annotation
  
    # read files (fread)
            # option -- multiple files or ONE annotated file (one annotated file requires AT LEAST a column indicating difference between samples)
          
    # merge
    # add groups
    # add keywords
    # remove unwanted columns
    # subsample
    # remove duplicates

    # add cell counts per sample

  # pre-processing

    # compensation?
    # transformation
    # batch-correction

        ## --> important interative step


  # Analysis prep
    # specify valid cellular parameters
    # specify stable phenotypers
    # specify downsampling targets AND METHOD prior to dim. red.

  ########## END USER INPUT ############

  # Clustering

    # Run FlowSOM n x times
    # Run something else n x times


  # Generate summary data, heatmaps, and autographs

    # Create and save sum.tables
    # Generate heatmap.cellnum
    # Generate heatmap.cellnum (fold)
    # Generate heatmap.mfi
    # Generate heatmap.mfi(fold)

    # Generate AutoGraphs (by population)
    # Generate AutoGraphs (by MFI)


  # Dim red

    # Downsample (using specified method)
    # Run tSNE -- add tSNE X and Y -- cycle, tSNE2.x and tSNE2.y etc
        # Print tSNEplots (option in function, yes/no)
    # Run UMAP
        # Print UMAP plots (option in function, yes/no)


  # Dim red plotting

    # tSNEplots -- cycle through combinations of tSNE/UMAPs with X and Y

      # if numeric, plot with color
      # if factor, plot with other colour scheme and LABEL



  #### THEN GO AWAY AND ASSESS

  # Read data (single file)

  # cluster.annotate
  
  # sumtables

  # heatmaps

  # autograph

  # tSNEplots





### CytoTools::read.files
    
    #' Read files function
    #'
    #' This function allows you to read in sample files (.csv) into a list
    #' @param file.loc What is the location of your files? ?? efaults to TRUE ??
    #' @keywords cats
    #' @export
    #' @examples
    #' cat_function()

CytoTools.read.files <- function(file.loc){
    
    file.names <- list.files(path=file.loc, pattern = ".csv")
    
    as.matrix(file.names)
    
      ## Read data from Files into list of data frames
      DataList=list() # Creates and empty list to start 
      Length_check = list() # creates an empty list
      ColName_check = list() 
      nrow_check = list()
    
    for (file in file.names) { # Loop to read files into the list
      tempdata <- fread(file, check.names = FALSE)
      file <- gsub(".csv", "", file)
      DataList[[file]] <- tempdata
      }
    
    for(i in c(1:(length(DataList)))){Length_check[[i]] <- length(names(DataList[[i]]))} # creates a list of the number of columns in each sample
    for(i in c(1:(length(DataList)))){ColName_check[[i]] <- names(DataList[[i]])}
    name.table <- data.frame(matrix(unlist(ColName_check), nrow = length(DataList), byrow = T))
    for(i in c(1:(length(DataList)))){nrow_check[[i]] <- nrow(DataList[[i]])}
    
    rm(tempdata)
    
    assign("DataList", DataList, envir = globalenv())
    assign("Length_check", Length_check, envir = globalenv())
    assign("ColName_check", ColName_check, envir = globalenv())
    assign("nrow_check", nrow_check, envir = globalenv())
    
    #return(DataList)
    msg <- "Files have been imported into a list called 'DataList"
    print(msg)
    
    }


### CytoTools::merge.csvs
    
CytoTools.merge.files <- function(x){
      
      ## Create a list of 'SampleName' and SampleNo' entries -- 1:n of samples, these will be matched in order
      AllFileNames <- c(names(x))
      AllFileNames # Character (words)
      
      AllFileNos <- c(1:(length(x)))       ### <-- changed to 4+ to match with sample
      AllFileNos # Intiger/numeric (numbers)
      
      ## Add 'SampleNo' and SampleName' to each 
      for (i in AllFileNos) {x[[i]]$SampleNo <- i}
      for (a in AllFileNames) {x[[a]]$SampleName <- a} # Inserting names doesn't stop the rest working, just messes with some auto visual stuff
      
      ###
      
      cell.dat <- rbind.fill(x) 
      cell.dat
      
      assign("cell.dat", cell.dat, envir = globalenv())
      assign("AllFileNames", AllFileNames, envir = globalenv())
      assign("AllFileNos", AllFileNos, envir = globalenv())
      
      return(head(cell.dat))
    }


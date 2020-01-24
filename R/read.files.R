#' read.files - Function to read data from CSV or FCS files into a list
#'
#' This function allows you to read in sample files (.csv or .fcs) into a list called 'data.list'.
#'
#' @param file.loc What is the location of your files? Defaults to current working directory (getwd()).
#' @param file.type What type of files do you want to read. Can be ".csv" or ".fcs". Defaults to ".csv".
#' @param do.embed.file.names Do you want to embed each row (cell) of each file with the name name? Defaults to TRUE.
#'
#' @return a list (called data.list) of dataframes -- one file per dataframe.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
#'
#' @examples
#' read.files()
#'
#' @export

    read.files <- function(file.loc = getwd(),
                           file.type = ".csv",
                           do.embed.file.names = TRUE){

        ## testing
            #file.loc <- getwd()
            #file.type <- ".csv"
            #do.embed.file.names = TRUE

        ## Read data from Files into list of data frames
        data.list=list() # Creates and empty list to start
        ncol.check = list() # creates an empty list
        colName.check = list()
        nrow.check = list()

        ## For reading CSV files
        if(file.type == ".csv"){
          file.names <- list.files(path=file.loc, pattern = file.type)

          for (file in file.names) { # Loop to read files into the list
            tempdata <- fread(file, check.names = FALSE)
            file <- gsub(".csv", "", file)
            data.list[[file]] <- tempdata
          }

          rm(tempdata)
          msg <- "CSV files have been imported into a list called 'data.list"
        }

        ## For reading FCS files
        if(file.type == ".fcs"){

          file.names <- list.files(path=file.loc, pattern = file.type)

          for (file in file.names) { # Loop to read files into the list
            tempdata <- exprs(read.FCS(file, transformation = FALSE))
            tempdata <- tempdata[1:nrow(tempdata),1:ncol(tempdata)]
            tempdata <- as.data.table(tempdata)
            file <- gsub(".fcs", "", file)
            data.list[[file]] <- tempdata
            #data.list[[file]] <- as.data.frame(data.list[[file]])
          }

          rm(tempdata)
          msg <- "FCS files have been imported into a list called 'data.list"
        }


        if(do.embed.file.names == TRUE){

          ## Create a list of 'SampleName' and SampleNo' entries -- 1:n of samples, these will be matched in order
          all.file.names <- c(names(data.list))
          all.file.names # Character (words)

          all.file.nums <- c(1:(length(data.list)))       ### <-- changed to 4+ to match with sample
          all.file.nums # Intiger/numeric (numbers)

          ## Add 'SampleNo' and SampleName' to each
          for (a in all.file.names) {data.list[[a]]$FileName <- a} # Inserting names doesn't stop the rest working, just messes with some auto visual stuff
          for (i in all.file.nums) {data.list[[i]]$FileNo <- i}

          assign("all.file.names", all.file.names, envir = globalenv())
          assign("all.file.nums", all.file.nums, envir = globalenv())
        }


        for(i in c(1:(length(data.list)))){ncol.check[[i]] <- length(names(data.list[[i]]))} # creates a list of the number of columns in each sample
        for(i in c(1:(length(data.list)))){colName.check[[i]] <- names(data.list[[i]])}
        name.table <- data.frame(matrix(unlist(colName.check), nrow = length(data.list), byrow = T))
        for(i in c(1:(length(data.list)))){nrow.check[[i]] <- nrow(data.list[[i]])}

        ## Save relevant values
        assign("name.table", name.table, envir = globalenv())

        ncol.check <- as.matrix(ncol.check)
        nrow.check <- as.matrix(nrow.check)

        assign("ncol.check", ncol.check, envir = globalenv())
        #assign("colName.check", colName.check, envir = globalenv())
        assign("nrow.check", nrow.check, envir = globalenv())

        #assign("data.list", data.list, envir = globalenv())
        return(data.list)

        ## Result message
        print(msg)
        }

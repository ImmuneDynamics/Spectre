#' read.files - Function to read data from CSV or FCS files into a list.
#'
#' This function allows you to read in sample files (.csv or .fcs) into a list, where each file is saved as a data.table.
#'
#' @usage read.files(file.loc, file.type, do.embed.file.names, header)
#'
#' @param file.loc DEFAULT = getwd(). What is the location of your files?
#' @param file.type DEFAULT = ".csv". What type of files do you want to read. Can be ".csv" or ".fcs".
#' @param do.embed.file.names DEFAULT = TRUE. Do you want to embed each row (cell) of each file with the name name?
#' @param header DEFAULT = TRUE. Does the first line of data contain column names?
#'
#' @return Returns a list of data.tables -- one per CSV file.
#'
#' @author
#' Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#'
#' @references Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
#'
#' @examples
#' data.list <- read.files(file.loc = getwd(), file.type = ".csv", do.embed.file.names = TRUE)
#'
#' @import data.table
#'
#' @export

read.files <- function(file.loc = getwd(),
                       file.type = ".csv",
                       do.embed.file.names = TRUE,
                       header = TRUE)
{

    ## Check that necessary packages are installed
        if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
        if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')

    ## Require packages
        require(Spectre)
        require(data.table)

    ## Initial warnings
        orig_wd <- getwd()
        
        if (!dir.exists(paste(orig_wd, file.loc, sep='/')) &
            !dir.exists(file.loc)) {
            warning("We were not able to find the directory specified by file.loc. Are you sure that location exists?")
        }
        
        setwd(file.loc)
        wd <- getwd()

        if(length(list.files(path=wd, pattern = file.type)) == 0){
          warning("We did not find any files in that directory, are you sure this is the right place?")
        }

    ## Read data from Files into list of data frames
        data.list=list() # Creates and empty list to start

        ncol.check = list() # creates an empty list
        colName.check = list()
        nrow.check = list()

    ## For reading CSV files
        if(file.type == ".csv"){
          file.names <- list.files(path=wd, pattern = file.type)

          for (file in file.names) { # Loop to read files into the list
            tempdata <- data.table::fread(file, check.names = FALSE, header = header)
            file <- gsub(".csv", "", file)
            data.list[[file]] <- tempdata
          }

          rm(tempdata)
          msg <- "CSV files have been imported into a list"
        }

    ## For reading FCS files
        if(file.type == ".fcs"){
          if(!is.element('flowCore', installed.packages()[,1])) stop('flowCore is required but not installed')
          require(flowCore)
          
          file.names <- list.files(path=wd, pattern = file.type)

          for (file in file.names) { # Loop to read files into the list
            tempdata <- exprs(flowCore::read.FCS(file, transformation = FALSE))
            tempdata <- tempdata[1:nrow(tempdata),1:ncol(tempdata)]
            tempdata <- as.data.table(tempdata)
            file <- gsub(".fcs", "", file)
            data.list[[file]] <- tempdata
            #data.list[[file]] <- as.data.frame(data.list[[file]])
          }

          rm(tempdata)
          msg <- "FCS files have been imported into a list"
        }

    ## For embedding file names

        if(do.embed.file.names == TRUE){

          ## Create a list of 'SampleName' and SampleNo' entries -- 1:n of samples, these will be matched in order
          all.file.names <- c(names(data.list))
          all.file.names # Character (words)

          all.file.nums <- c(1:(length(data.list)))       ### <-- changed to 4+ to match with sample
          all.file.nums # Intiger/numeric (numbers)

          ## Add 'SampleNo' and SampleName' to each
          for (a in all.file.names) {data.list[[a]]$FileName <- a} # Inserting names doesn't stop the rest working, just messes with some auto visual stuff
          for (i in all.file.nums) {data.list[[i]]$FileNo <- i}

        }

    ## Return result and result message
        setwd(orig_wd)
        
        return(data.list)
        message(msg)
}

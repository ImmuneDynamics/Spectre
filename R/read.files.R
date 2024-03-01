#' read.files
#' 
#' Function to read cytometry data stored in either CSV or FCS files into 
#' either a list or a Spectre object.
#' 
#' See details for more information.
#'
#' @details
#' If you choose the output as a list, the function will read each file
#' into a data.table and store it in an element of the list.
#' Essentially 1 element = 1 file.
#' 
#' If you choose the output as Spectre object, the `cytometry_data` slot
#' will be filled with a data.table containing the concatenation of all the files.
#' A new column will be created to denote the file the cell comes from.
#'
#' @usage read.files(file.loc)
#'
#' @param file.loc What is the location of your files?
#' @param file.type DEFAULT = ".csv". What type of files do you want to read.
#' Can be ".csv" or ".fcs".
#' @param nrows DEFAULT = NULL. Can specify a numerical target for the number of
#' cells (rows) to be read from each file. 
#' Please note, order is random in FCS files.
#' @param do.embed.file.names DEFAULT = TRUE. 
#' Do you want to embed each row (cell) of each file with the name name?
#' Only used if output is a list. 
#' @param header DEFAULT = TRUE. Does the first line of data contain column names?
#' Only used if file.type is .csv.
#' @param truncate_max_range DEFAULT = TRUE. Whether to truncate the extreme
#' positive value to the instrument measurement range. 
#' Only used when reading from FCS files.
#' @param as_spectre_object DEFAULT = FALSE. 
#' Whether to return the files as a Spectre object.
#' @param verbose DEFAULT = FALSE.
#' If TRUE, the function will print progress updates as it executes.
#'
#' @return Either a list of data.tables (one element per file) or a SpectreObject
#'
#' @author Thomas M Ashhurst \email{thomas.ashhurst@@sydney.edu.au},
#' Felix Marsh-Wakefield \email{felix.marsh-wakefield@@sydney.edu.au},
#' Givanna Putri
#'
#' @references Ashhurst, T. M., et al. (2019).
#'   \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
#'
#' @examples
#' data.list <- read.files(file.loc = getwd(), file.type = ".csv", do.embed.file.names = TRUE)
#'
#'
#' @export
read.files <- function(file.loc,
                       file.type = c(".csv", ".fcs"),
                       nrows = NULL,
                       do.embed.file.names = TRUE,
                       truncate_max_range = TRUE,
                       header = TRUE,
                       as_spectre_object = FALSE,
                       verbose = FALSE) {
    if (!dir.exists(file.loc)) {
        error(paste(file.loc, "directory does not exist!"))
    }
    
    file.type <- match.arg(file.type)
    file.names <- list.files(
        path = file.loc, 
        pattern = paste0("\\", file.type, "$")
    )
    
    if (length(file.names) == 0) {
        error("We did not find any files in that directory, are you sure this is the right place?")
    }
    
    
    ## For reading CSV files
    
    if (file.type == ".csv") {
        
        data.list <- lapply(seq(1, length(file.names)), function(file_no) {
            file <- file.path(file.loc, file.names[[file_no]])
            
            if (is.null(nrows)) {
                
                if (verbose) {
                    message(paste("Reading file", file))
                }
                
                tempdata <- data.table::fread(file, check.names = FALSE, header = header)
            } else {
                
                if (verbose) {
                    message(paste("Reading only the first", nrows, 
                                  "rows (cells) for file",
                                  file))
                }
                
                tempdata <- data.table::fread(file,
                                              check.names = FALSE,
                                              header = header,
                                              nrows = nrows
                )
            }
            
            if (do.embed.file.names) {
                tempdata$FileName <- gsub("\\.csv$", "", file.names[[file_no]])
                tempdata$FileNo <- file_no
            }
            
            
            return(tempdata)
        })
        
        names(data.list) <- gsub("\\.csv$", "", file.names)
        
        msg <- "CSV files have been imported into a list"
    }
    
    
    ## For reading FCS files
    
    else if (file.type == ".fcs") {
        
        data.list <- lapply(seq(1, length(file.names)), function(file_no) {
            file <- file.path(file.loc, file.names[[file_no]])
            
            if (is.null(nrows)) { 
                ## If nrows not specified
                if (verbose) {
                    message(paste("Reading file", file))
                }
                
                x <- flowCore::read.FCS(file,
                                        transformation = FALSE,
                                        truncate_max_range = truncate_max_range
                )
            }
            
            else { 
                ## If nrows specified
                if (verbose) {
                    message(paste("Reading only the first", nrows, 
                                  "rows (cells) for file",
                                  file))
                }
                x <- flowCore::read.FCS(file,
                                        transformation = FALSE,
                                        which.lines = nrows,
                                        truncate_max_range = truncate_max_range
                )
            }
            
            if (verbose) 
                message("Extracting data from FlowCore object")
            nms <- vector()
            for (o in c(1:nrow(x@parameters@data))) {
                pr <- x@parameters@data$name[[o]]
                st <- x@parameters@data$desc[[o]]
                
                if (!is.na(st)) {
                    nms <- c(nms, paste0(pr, "_", st))
                } else {
                    nms <- c(nms, pr)
                }
            }
            
            if (verbose) 
                message("Creating data.table object")
            tempdata <- exprs(x)
            tempdata <- tempdata[1:nrow(tempdata), 1:ncol(tempdata)]
            tempdata <- data.table::as.data.table(tempdata)
            names(tempdata) <- nms
            
            if (do.embed.file.names) {
                tempdata$FileName <- gsub("\\.fcs$", "", file.names[[file_no]])
                tempdata$FileNo <- file_no
            }
            
            return(tempdata)
        })
        
        names(data.list) <- gsub("\\.fcs$", "", file.names)
        
        msg <- "FCS files have been imported into a list"
    }
    
    if (verbose)
        message(msg)
    
    if (as_spectre_object) {
        
        if (verbose)
            message("Converting list to a Spectre object")
        
        data.list <- do.merge.files(data.list)
        obj <- SpectreObject(cytometry_data=data.list, citeseq_data=data.table())
        return(obj)
    } else {
        return(data.list)
    }
}
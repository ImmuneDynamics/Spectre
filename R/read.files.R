#' Read CSV or FCS files into a list.
#'
#' This function allows you to read in sample files (.csv or .fcs) into a list, 
#' where each file is saved as a data.table object.
#'
#' @param file.loc DEFAULT = getwd(). What is the location of your files?
#' @param file.type DEFAULT = ".csv". What type of files do you want to read. Can be ".csv" or ".fcs".
#' @param files DEFAULT = NULL. A vector of selected file names to import.
#' @param nrows DEFAULT = NULL. Can specify a numerical target for the number of cells (rows) to be read from each file. Please note, order is random in FCS files.
#' @param do.embed.file.names DEFAULT = TRUE. Do you want to embed each row (cell) of each file with the name name?
#' @param header DEFAULT = TRUE. Does the first line of data contain column names?
#' Only works for file.type = ".csv".
#'
#' @return Returns a list of data.tables -- one per file.
#'
#' @author
#' Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#' Givanna Putri
#'
#' @references Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
#'
#' @examples
#' # download sample data
#' download.file(url='https://github.com/ImmuneDynamics/data/blob/main/msCNS.zip?raw=TRUE', destfile = 'msCNS.zip', mode = 'wb')
#' unzip(zipfile = 'msCNS.zip')
#' data.list <- read.files(file.loc = "msCNS/data", file.type = ".csv", do.embed.file.names = TRUE)
#'
#' @import data.table
#' @import fs
#' @import flowCore
#'
#' @export
#' 
read.files <- function(
    file.loc = getwd(),
    file.type = c(".csv", ".fcs"),
    files = NULL,
    nrows = NULL,
    do.embed.file.names = TRUE,
    header = TRUE
) {
    file.type <- tryCatch(
        match.arg(file.type),
        error = function(e) {
            stop("Invalid value for 'file.type': must be either '.csv' or '.fcs'.", call. = FALSE)
        }
    )
    
    # Convert to absolute path (resolves relative paths)
    abs.file.loc <- fs::path_abs(file.loc)
    
    # Check if it's an existing directory
    if (!fs::dir_exists(abs.file.loc)) {
        stop(
            paste0("Directory not found: '", file.loc, "'\nAre you sure that location exists?")
        )
    }

    # Check if there are any files in it?
    if (length(list.files(path = abs.file.loc, pattern = file.type)) == 0){
        stop(
            paste0("We did not find any ", file.type,  " files in ", file.loc, ".\nAre you sure this is the right place?")
        )
    }

    if (is.null(files)) {
        file.names <- list.files(path = abs.file.loc, pattern = file.type)
    } else {
        file.names <- files
    }

    if (file.type == ".csv") {
        message("Reading CSV files...")
        
        data.list <- lapply(seq_len(length(file.names)), function(i) {
            file.name <- file.names[i]
            message(paste("Reading", file.name))
            if (is.null(nrows)) {
                tempdata <- data.table::fread(
                    file.path(abs.file.loc, file.name),
                    check.names = FALSE,
                    header = header
                )
            } else {
                message(paste0("Reading ", nrows, " rows (cells) per file"))
                tempdata <- data.table::fread(
                    file.path(abs.file.loc, file.name),
                    check.names = FALSE,
                    header = header,
                    nrows = nrows
                )
            }
            if (do.embed.file.names) {
                tempdata$FileName <- gsub(".csv", "", file.name)
                tempdata$FileNo <- i
            }
            return(tempdata)
        })
        names(data.list) <- gsub(".csv", "", file.names)

    } else if (file.type == ".fcs") {
        data.list <- lapply(seq_len(length(file.names)), function(i) {
            file.name <- file.names[i]
            message(paste("Reading", file.name))
            if (is.null(nrows)) {
                x <- flowCore::read.FCS(
                    file.path(abs.file.loc, file.name), 
                    transformation = FALSE
                )
                
            } else {
                message(paste0("Reading ", nrows, " rows (cells) per file"))
                x <- flowCore::read.FCS(
                    file.path(abs.file.loc, file.name), 
                    transformation = FALSE, 
                    which.lines = nrows
                )
            }
            channel_name <- x@parameters@data$name
            antibody_name <- x@parameters@data$desc
            tempdata_colnames <- ifelse(
                is.na(antibody_name),
                channel_name,
                paste0(channel_name, "_", antibody_name)
            )
            
            tempdata <- as.data.table(exprs(x))
            names(tempdata) <- tempdata_colnames
            
            if (do.embed.file.names) {
                tempdata$FileName <- gsub(".fcs", "", file.name)
                tempdata$FileNo <- i
            }

            return(tempdata)
        })
        names(data.list) <- gsub(".fcs", "", file.names)
    }
    return(data.list)
}

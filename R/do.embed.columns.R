#' do.embed.columns
#'
#' @usage do.embed.columns(x, type, base.name, col.name, match.to, new.cols, rmv.ext)
#'
#' @param x NO DEFAULT. A list of dataframes (samples) that will recieve the new columns.
#' @param type DEFAULTS TO "data.table". Can be "data.table" or "list". Whether you wish to embed new values into a single dataframe based on values in another column, or dataframes as part of a list based on the name of the dataframes.
#' @param base.name NO DEFAULT. Name or number of the column containing values that new values will be matched to. Only required when embedding into a single data.table/data.frame -- when embedding into a list of data.tables/data.frames, matching is done based on the name of each data.table/data.frame
#' @param col.name NO DEFAULT. Character, desired name for the new column.
#' @param match.to NO DEFAULT. Avector that specifies which terms are used for matching (e.g. c("Sample 1", "Sample 2", "Sample 3", "Sample 4"), or c("Cluster 1", "Cluster 2"), etc).
#' @param new.cils NO DEFAULT. A vector of new values to embed as a new column, matched to the corresponding value in 'match.to' (e.g. c("Group 1", "Group 1", "Group 2", "Group 2"), or c("T cells", "B cells"), etc).
#' @param rmv.ext DEFAULTS TO TRUE. Logical, can be TRUE or FALSE. Removes the ".csv" or ".fcs" extension from a the 'match.to' vector -- especially useful if 'match.to' is a list of sample names that end in .csv or .fcs.
#'
#' @return For type = "data.table", returns the data.table with new columns embedded. For type = "list", returns list with newly embedded values in each data.table within the list.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @examples
#' # For type = "list"
#' data.list <- .........
#'
#' # For type = "data.table"
#' cell.dat <- ........
#'
#' @export

embed.columns <- function(x, # the list of dataframes (samples) where each dataframe will have the columns embedded
                          type = "data.table", # "list" or "data.frame"
                          base.name, #
                          col.name,
                          match.to, # vector of the element names to match to
                          new.cols, # vector of the new values to embedded
                          rmv.ext = TRUE
                          )
{

    ##############################################
    ### For embedding into data.tables in a list
    ##############################################

    ## Test for list

        # x <- data.list
        # type <- "list"
        # col.name <- names(meta.dat$sampleDetails[c(2)])
        # match.to <- meta.dat$sampleDetails[c(1)]
        # new.cols <- meta.dat$sampleDetails[c(2)]
        # rmv.ext <- FALSE # will remove the .fcs or .csv file extension in the sample details df, if it is present.

    if(type == "list"){

      if(isFALSE(length(match.to) == length(new.cols))){  # if lengh doesn't match
        stop("Error - length of match.to and new.cols does not match.")
      }

      if(length(match.to) == length(new.cols)){

        ##
        cols <- cbind(match.to, new.cols)
        cols <- as.data.frame(cols)

        cols[1] <- lapply(cols[1], function(cols){cols <- gsub(".csv", "", cols); cols})
        cols[1] <- lapply(cols[1], function(cols){cols <- gsub(".fcs", "", cols); cols})

        ## Check
        chk <- any(isFALSE(cols[1] == names(x)))

        ##
        if(chk == FALSE){
          for(i in c(1:nrow(cols[1]))){
            nme <- cols[i,1]
            x[[nme]][[col.name]] <- NA
            x[[nme]][[col.name]] <- cols[i,2]
          }
        }

        if(chk == TRUE){
          stop("Error - list element list and 'match.to' do not match")
        }
      }
      #assign("data.list", x, envir = globalenv())
      return(x)
    }

  ##############################################
  ### For single data.tables
  ##############################################

  ## Test for data.table

      # x <- demo.clustered
      # type <- "data.table"
      # base.name <- "FlowSOM_metacluster"
      # col.name <- "PopName"
      #
      #     unique(x$FlowSOM_metacluster)
      #     min(x$FlowSOM_metacluster)
      #     max(x$FlowSOM_metacluster)
      #
      # match.to <- c(1:40)
      # new.cols <- c(rep("PMN", 20), rep("T cell", 20))
      # rmv.ext <- FALSE # will remove the .fcs or .csv file extension in the sample details df, if it is present.

  if(type == "data.table"){

    if(isFALSE(length(match.to) == length(new.cols))){  # if lengh doesn't match
      stop("Error - length of match.to and new.cols does not match.")
    }

    if(length(match.to) == length(new.cols)){

      ##
      cols <- cbind(match.to, new.cols)
      cols <- as.data.frame(cols)

      cols[1] <- lapply(cols[1], function(cols){cols <- gsub(".csv", "", cols); cols})
      cols[1] <- lapply(cols[1], function(cols){cols <- gsub(".fcs", "", cols); cols})

      ## Check
      chk <- any(isFALSE(cols[1] == names(x)))

      ##
      if(chk == FALSE){
        x[[col.name]] <- NA
        for(i in (cols[,1])){
          x[[col.name]][x[[base.name]] == i] <- as.vector(cols[i,2])
        }
          embed.res <- as.data.frame(x[[col.name]])
          names(embed.res) <- col.name
      }

      if(chk == TRUE){
        stop("Error - list element list and 'match.to' do not match")
      }
    }
    #assign("embed.res", embed.res, envir = globalenv())
    #return(embed.res)
    return(x)
  }

}


#' embed.columns
#'
#' @usage embed.columns(x, ...)
#'
#' @param x NO DEFAULT. A list of dataframes (samples) that will recieve the new columns.
#' @param type DEFAULTS TO "data.frame". Can be "data.frame" or "list". Whether you wish to embed new values into a single dataframe based on values in another column, or dataframes as part of a list based on the name of the dataframes.
#' @param base.name NO DEFAULT. The column containing values that new values will be matched to. Only required when embedding into a single data frame -- when embedding into a list of dataframes, matching is done based on the name of each dataframe.
#' @param col.name NO DEFAULT. Character, desired name for the new column.
#' @param match.to NO DEFAULT. A character vector that specifies which dataframe in the list is used for matching (e.g. "Sample 1", "Sample 2", "Sample 3", "Sample 4", etc).
#' @param new.cils NO DEFAULT. A vector of new values to embed to each dataframe, matched to the corresponding value in 'match.to' (e.g. "Group 1", "Group 1", "Group 2", "Group 2", etc).
#' @param rmv.ext DEFAULTS TO TRUE. Logical, remove the ".csv" or ".fcs" extension from a the 'match.to' vector -- especially useful if 'match.to' is a list of sample names that end in .csv or .fcs.
#'
#' @return For type = "data.frame", returns a vector of the new row values (called embed.res), that can be attached to the starting dataset using cbind. For type = "list", embeds new columns into each row (cell) of dataframes in a list called 'data.list'.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @examples embed.columns()
#'
#' @export

embed.columns <- function(x, # the list of dataframes (samples) where each dataframe will have the columns embedded
                          type, # "list" or "data.frame
                          base.name, #
                          col.name,
                          match.to, # vector of the element names to match to
                          new.cols, # vector of the new values to embedded
                          rmv.ext = TRUE
                          )
{

    ### Test data

     ## Test for list

        # x <- data.list
        # type <- "list"
        # col.name <- names(meta.dat$sampleDetails[c(2)])
        # match.to <- meta.dat$sampleDetails[c(1)]
        # new.cols <- meta.dat$sampleDetails[c(2)]
        # rmv.ext <- FALSE # will remove the .fcs or .csv file extension in the sample details df, if it is present.

      ## Test for data frame

        # x <- cell.dat
        # type <- "data.frame"
        # base.name <- "FlowSOM_metacluster"
        # col.name <- "PopName"
        #
        #     unique(cell.dat$FlowSOM_metacluster)
        #     min(cell.dat$FlowSOM_metacluster)
        #     max(cell.dat$FlowSOM_metacluster)
        #
        # match.to <- c(1:40)
        # new.cols <- c(rep("PMN", 20), rep("T cell", 20))
        # rmv.ext <- FALSE # will remove the .fcs or .csv file extension in the sample details df, if it is present.


    ##############################################
    ### For embedding into dataframes in a list
    ##############################################

    if(type == "list"){

      if(isFALSE(length(match.to) == length(new.cols))){  # if lengh doesn't match
        print("Error - length of match.to and new.cols does not match.")
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
          print("Error - list element list and 'match.to' do not match")
        }
      }
      #assign("data.list", x, envir = globalenv())
      return(data.list)
    }

  ##############################################
  ### For single data.frames
  ##############################################

  if(type == "data.frame"){

    if(isFALSE(length(match.to) == length(new.cols))){  # if lengh doesn't match
      print("Error - length of match.to and new.cols does not match.")
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
        print("Error - list element list and 'match.to' do not match")
      }
    }
    #assign("embed.res", embed.res, envir = globalenv())
    return(embed.res)
  }

}


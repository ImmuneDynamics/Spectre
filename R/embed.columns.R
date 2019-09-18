#' embed.columns
#'
#' @usage embed.columns(x, ...)
#'
#' @param x NO DEFAULT. A list of dataframes (samples) that will recieve the new columns.
#' @param type DEFAULTS TO "list". Can be "list", "data.frame" will be available in next release. Whether you wish to embed new values into dataframes as part of a list based on the name of the dataframes, or into a single dataframe based on values in another column.
#' @param col.name NO DEFAULT. Character, desired name for the new column.
#' @param match.to NO DEFAULT. A character vector that specifies which dataframe in the list is used for matching (e.g. "Sample 1", "Sample 2", "Sample 3", "Sample 4", etc).
#' @param new.cils NO DEFAULT. A vector of new values to embed to each dataframe, matched to the corresponding value in 'match.to' (e.g. "Group 1", "Group 1", "Group 2", "Group 2", etc).
#' @param rmv.ext DEFAULTS TO TRUE. Logical, remove the ".csv" or ".fcs" extension from a the 'match.to' vector -- especially useful if 'match.to' is a list of sample names that end in .csv or .fcs.
#'
#' @return Embeds new columns into each row (cell) of dataframes in a list
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
                          col.name,
                          match.to, # vector of the element names to match to
                          new.cols, # vector of the new values to embedded
                          rmv.ext = TRUE
                          )
{

    ### Test data

        # x <- data.list
        # type <- "list"
        # col.name <- names(meta.dat$sampleDetails[c(2)])
        # match.to <- meta.dat$sampleDetails[c(1)]
        # new.cols <- meta.dat$sampleDetails[c(2)]
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
      assign("data.list", x, envir = globalenv())
    }

  ##############################################
  ### For single data.frames
  ##############################################

  if(type == "data.frame"){
    print("embedding into a single data.frame not supported in this version")
  }

}




# if(type == "data.frame"){
#   print("Embedding columns into a single dataframe is not yet supported")
# }


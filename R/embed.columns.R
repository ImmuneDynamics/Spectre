#' embed.columns
#'
#' @usage embed.columns(x, ...)
#'
#' @return Embeds new columns into samples
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}. Helpful examples at \url{https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/}
#'
#' @examples embed.columns()
#'
#' @export

embed.columns(x, # the list of dataframes (samples) where each dataframe will have the columns embedded
              new.cols, # a table of filenames and new columns to add
              file.col  # which column contains the reference filenames
              )
{

  ### Test data

      # x <- data.list
      # new.cols <- sample.table[c(1:4)]
      # file.col <- "Filename"

  ### Sorting out the columns

      # Remove any '.fcs'
      fls <- as.matrix(new.cols[file.col])
      for(i in c(1:length(fls))){
        test <- fls[i,]
        test
        test <- gsub(".fcs", "", test)
        fls[i,] <- test
      }

      # Remove any '.csv'
      fls <- as.matrix(new.cols[file.col])
      for(i in c(1:length(fls))){
        test <- fls[i,]
        test
        test <- gsub(".csv", "", test)
        fls[i,] <- test
      }
      fls

  ### Remove filename from 'new.cols'
      new.cols <- new.cols[-(which( colnames(new.cols)==file.col ))]

  ### Check filenames are consistent

      checks <- fls == names(x)
      all(checks)

  ### Perform embedding

      if(all(checks) == TRUE){ ## Embed metadata
        for(i in c(1:length(names(x)))){
          #i <- 1
          for(a in names(new.cols)){
            x[[i]][[a]] <- NA # fills a new colum
            x[[i]][[a]] <- new.cols[i,a]
          }
        }
      }

      if(all(checks) == FALSE){ ## Warning message
        print("The list of file names in the list of samples, and the file names on your table to do not match. Unable to embed new columns"
        )
      }

}

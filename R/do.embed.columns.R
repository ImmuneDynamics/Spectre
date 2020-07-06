#' do.embed.columns
#'
#' @usage do.embed.columns(dat, base.col, add.dat, add.by, rmv.ext)
#'
#' @param dat NO DEFAULT. A data.table (or data.frame) containing the data to have new values added to.
#' @param base.col NO DEFAULT. Name of the column containing values that new values will be matched to.
#' @param add.dat NO DEFAULT. A data table of new values to embed as a new columns, with one column containing values used for matching with the target data.table.
#' @param add.by NO DEFAULT. Character, name of the column in add.dat that is used for matching to the taret dataset.
#' @param rmv.ext DEFAULTS TO TRUE. Logical, can be TRUE or FALSE. Removes the ".csv" or ".fcs" extension from a the 'match.to' vector -- especially useful if 'match.to' is a list of sample names that end in .csv or .fcs.
#'
#' @return Returns the data.table with new columns embedded..
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @examples
#'
#' @export

do.embed.columns <- function(dat, # the list of dataframes (samples) where each dataframe will have the columns embedded
                             base.col, # column name in the actual dataset

                             add.dat,
                             add.by,
                             rmv.ext = TRUE)
{
  dat <- as.data.table(dat)
  add.dat <- as.data.table(add.dat)

  if(rmv.ext == TRUE){
    for(i in c(1:length(names(add.dat)))){
      if(is.numeric(add.dat[[i]]) == FALSE){
        temp <- add.dat[[i]]
        temp <- gsub("*.csv", "", temp)
        temp <- gsub("*.fcs", "", temp)
        add.dat[[i]] <- temp
      }
    }
  }

  add.dat <- cbind(add.dat[,add.by,with = FALSE], add.dat[,!add.by,with = FALSE])
  names(add.dat)[c(1)] <- base.col

  if(class(dat[[base.col]]) != class(add.dat[[base.col]])){
    warning(paste0("The column '", base.col, "' in your main dataset", " is ", class(dat[[base.col]]), ", whereas the column '", add.by, "' in the 'to add' data is ", class(add.dat[[base.col]]), ". You may need to adjust their types to ensure they are the same. If you are matching based on character values (filenames, popualtion names etc, then both need to be character."))
  }

  res = merge(dat, add.dat, by = base.col, sort = FALSE)

  return(res)
}


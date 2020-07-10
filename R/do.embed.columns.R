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

  ## Packages
  if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
  if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')

  ### Require packages
  require(Spectre)
  require(data.table)

  ## Test data
      # dat <- Spectre::demo.start
      # base.col <- "FileName"
      #
      # add.dat <- data.table(Filename = unique(dat$FileName),
      #                       NewColA = c('A','B','C','D','E','F','G','H','I','J','K', 'L'),
      #                       NewColB = c(1,2,3,4,5,6,7,8,9,10,11,12))
      #
      # add.by <- "Filename"
      # rmv.ext = TRUE

  ##

  dat.names <- names(dat)
  add.dat.names <- names(add.dat)

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

  dat <- as.data.table(dat)
  add.dat <- as.data.table(add.dat)

  p1 <- add.dat[,add.by,with = FALSE]
  p2 <- add.dat[,setdiff(names(add.dat),add.by),with = FALSE]

  added.names <- names(p2)

  add.dat <- cbind(p1, p2)

  rm(p1)
  rm(p2)

  names(add.dat)[c(1)] <- base.col

  dat[[base.col]] <- as.factor(dat[[base.col]])
  add.dat[[base.col]] <- as.factor(add.dat[[base.col]])

  # class(dat[[base.col]])
  # class(add.dat[[base.col]])

  if(class(dat[[base.col]]) != class(add.dat[[base.col]])){
    warning(paste0("The column '", base.col, "' in your main dataset", " is ", class(dat[[base.col]]), ", whereas the column '", add.by, "' in the 'to add' data is ", class(add.dat[[base.col]]), ". You may need to adjust their types to ensure they are the same. If you are matching based on character values (filenames, popualtion names etc, then both need to be character."))
  }

  dat.names
  added.names

  res.1 = merge(dat, add.dat, by = base.col, sort = FALSE)
  res.1 <- res.1[,c(dat.names,added.names),]

  #res.2 <- res.1[,add.dat.names[c(2:length(add.dat.names))], with = FALSE]

  #res.3 <- cbind(dat, res.2)

  # qrs <- dat[,base.col,with = FALSE]
  #
  # setkeyv(qrs, base.col)
  # setkeyv(add.dat, base.col)
  #
  # qrs[,list(add.dat,qrs)]
  # add.dat[,list(qrs,add.dat)]
  #
  # setkey(dat, NULL)
  # setkeyv(add.dat, NULL)

  #return(res.3)
  return(res.1)
}


#' do.normalise - Normalise data in selected columns between two values, usually 0 and 1.
#'
#' This function allows you to normalise the data in selected columns between two values, usually 0 and 1.
#'
#' @seealso \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @param dat NO DEFAULT. data.table Input sample.
#' @param use.cols NO DEFAULT. Vector of character column names -- these columns will be normalised and added to the data.table as new columns.
#' @param new.min DEFAULT = 0. The new minimum value.
#' @param new.max DEFAULT = 1. The new maximum value.
#' @param zero.drops DEFAULT = NULL. This is an experimental calculation which should reduce noise from negative values. Use with caution. Supplying a vector of values for each column (that matches the order or use.cols) that reflects the negative 'cutoff' for eachc column. Values below this will be turned into this value, which will all then be converted to 0.
#'
#' @return A data.table with new columns added, that contain the normalised data.
#'
#' @import data.table
#'
#' @export

do.normalise <- function(dat,
                         use.cols,
                         new.min = 0,
                         new.max = 1,
                         zero.drops = NULL){ # null by default

  ### Packages
  require('data.table')
  
  ### Checks
  if(!is.null(zero.drops)){
    message("Zero drop transformation (as a part of the normalisation function) is experimental - please use with caution.")
    if(length(use.cols) != length(zero.drops)){
      warning("The number of your columns to modify the number of cutoff values. Please check your input.")
    }
  }

  ### Create normalisation function
  #norm.fun <- function(x) {(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) -min(x, na.rm=TRUE))}
  norm.fun <- function(x) {(x - min(x))/(max(x)-min(x)) * (new.max - new.min) + new.min}
  #norm.fun <- function(x){scales::rescale(value, to=c(new.min,new.max))}

  ### Establish dataset
  value <- dat[,use.cols,with = FALSE]

  ### Changes values < cutoff to the cutoff values
  if(!is.null(zero.drops)){
    for(i in c(1:length(use.cols))){
      # i <- 1
      a <- use.cols[[i]]
      b <- zero.drops[[i]]
      temp <- value[,a,with = FALSE]
      temp[temp[[a]] < b,] <- b
      value[,a] <- temp
    }
  }

  ### Normalise between 0 to 1
  res <- as.data.table(lapply(value, norm.fun)) # by default, removes the names of each row

  if(!is.null(zero.drops)){
    names(res) <- paste0(names(res), "_norm_zd")
    dat <- cbind(dat, res)
  }

  if(is.null(zero.drops)){
    names(res) <- paste0(names(res), "_norm")
    dat <- cbind(dat, res)
  }

  ###
  return(dat)
}


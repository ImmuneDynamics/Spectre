#' Noise reduction
#'
#' Method to reduce noise
#' 
#' @export

do.noise.reduce <- function(dat,
                            use.cols,
                            cutoffs){ # null by default
  
  ### Checks
  message("Limit adjustment transformation is experimental - please use with caution.")
  if(length(use.cols) != length(cutoffs)){
    warning("The number of your columns to modify the number of cutoff values. Please check your input.")
  }
  
  ### Establish dataset
  value <- dat[,use.cols,with = FALSE]
  
  ### Changes values < cutoff to the cutoff values
  if(!is.null(cutoffs)){
    for(i in c(1:length(use.cols))){
      # i <- 1
      a <- use.cols[[i]]
      b <- cutoffs[[i]]
      temp <- value[,a,with = FALSE]
      temp[temp[[a]] < b,] <- b
      value[,a] <- temp
    }
  }
  
  ### Names
  names(value) <- paste0(names(value), "_noiseRed")
  dat <- cbind(dat, value)
  
  ###
  return(dat)
}

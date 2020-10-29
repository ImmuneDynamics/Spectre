#' do.mrn - Multipurpose function for arcsinh transformation, Reduction of noise, and Normalisation
#'
#' This function facilitates multiple data transformation functions in a single process, including arcsinh transformation, noise reduction, and normalisation
#'
#' @usage do.mrn(dat, use.cols)
#'
#' @param dat NO DEFAULT. A data.table
#' @param use.cols NO DEFAULT. Vector of column names that you wish to transform
#' @param asinh.cf DEFAULT = 15. The co-factor to use for arcsinh transformation
#' @param low.threshold DEFAULT = 0. Values below this number will be converted to this number.
#' @param norm DEFAULT = c(0,1). New MIN and MAX values for each channel, as c(MIN, MAX). 
#' 
#' @return Returns a data.table with transformed data added as new columns.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
#'
#' @examples
#' cell.dat <- do.mrn()
#'
#' @export

do.mrn <- function(dat,
                   use.cols,
                   asinh.cf = 15,
                   low.threshold = 0,
                   norm = c(0,1)
){
  
  ### Extract data
  
      value <- dat[,use.cols,with = FALSE]
  
  ### Check for numeric values
      
      if(isFALSE(all(sapply(value, is.numeric)))){
        message("It appears that one column in your dataset is non numeric")
        print(sapply(value, is.numeric))
        stop("do.asinh stopped")
      }
  
  ### Arcsinh
  
      if(!is.null(asinh.cf)){
        value <- value / asinh.cf
        value <- asinh(value) # value <- log(value + sqrt(value^2 + 1))
        names(value) <- paste0(names(value), "_asinh")
        
        message("Arcsinh transformation complete")
      }

  ### Noise reduction
      
      if(!is.null(low.threshold)){
        for(a in names(value)){
          b <- low.threshold
          temp <- value[,a,with = FALSE]
          temp[temp[[a]] < b,] <- b
          value[,a] <- temp
        }
        
        names(value) <- paste0(names(value), "_NR")
        
        message("Noise reduction complete transformation complete")
      }
      
  ### Normalisation
  
      if(!is.null(norm)){
        norm.fun <- function(x) {(x - min(x))/(max(x)-min(x)) * (norm[2] - norm[1]) + norm[1]}
        value <- as.data.table(lapply(value, norm.fun)) # by default, removes the names of each row
        names(value) <- paste0(names(value), "_norm")
        
        message("Normalisation between ", norm[1], " and ", norm[2], " complete")
      }
      
  ### Wrap up
      
      value <- as.data.table(value)
      res <- cbind(dat, value)
      return(res)
  
}
  
 
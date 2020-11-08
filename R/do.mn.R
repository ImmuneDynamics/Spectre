#' do.mn - Multipurpose function for arcsinh transformation, quantile normalisation, noise reduction, and normalisation
#'
#' This function facilitates multiple data transformation functions in a single process, including arcsinh transformation, noise reduction, and normalisation
#'
#' @usage do.mn(dat, use.cols)
#'
#' @param dat NO DEFAULT. A data.table
#' @param use.cols NO DEFAULT. Vector of column names that you wish to transform
#' @param asinh.cf DEFAULT = 15. The co-factor to use for arcsinh transformation
#' @param norm DEFAULT = c(0,1). New MIN and MAX values for each channel, as c(MIN, MAX). 
#' @param lower.threshold DEFAULT = 0. Values below this number will be converted to this number.
#' @param upper.Q DEFAULT = 0.995. To not utilise the upper quantile function, simply use upper.Q = NULL or upper.Q = 1. When normalising, the new 'maximum' value (defined in the 'norm' argument, typically '1') will be set to this quantile target (for upper.Q = 0.995, this is the 99.5th percentile). As such, any values above the 99.5th percentile will increase above this new maximum value proportionally. Any target between 0 and 1 can be set. 
#' @param append.name DEFAULT = '_transf'. Text to be appended to resultant column name
#' @param info.dir DEFAULT = getwd(). Directory address to save a .txt. file with the settings used for the MRN function.
#' 
#' @return Returns a data.table with transformed data added as new columns.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
#'
#' @examples
#' cell.dat <- do.mn()
#'
#' @export

do.mn <- function(dat,
                  use.cols,
                   
                  asinh.cf = 15,
                  norm = c(0,1),
                  
                  lower.threshold = 0,
                  upper.Q = 0.995,
                  
                  append.name = '_transf',
                  save.info = TRUE,
                  info.dir = getwd()
){
  
  ### Libraries
  
      library(Spectre)
      library(data.table)
  
  ### Test data
      
      # dat <- Spectre::demo.start
      # use.cols <- names(dat)[c(2:10)]
      # 
      # asinh.cf <- 1000
      # norm = c(0,1)
      # 
      # lower.threshold = 0
      # upper.Q = 0.995
      # 
      # save.info = TRUE
      # append.name = '_transf'
      # info.dir = getwd()
  
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
        #names(value) <- paste0(names(value), "_asinh")
        
        message("Arcsinh transformation complete")
      }
 
  ### Lower threshold
      
      if(!is.null(lower.threshold)){
        for(a in use.cols){
          
          b <- lower.threshold
          temp <- value[,a,with = FALSE]
          temp[temp[[a]] < b,] <- b
          value[,a] <- temp
          
          rm(a)
          rm(b)
          rm(temp)
        }
        message("Lower value thresholding is complete")
      }
      
  ### Normalisation, including upper Q target
      
      if(!is.null(norm)){
        
        ## No upper Q (upper Q is NULL)
        if(is.null(upper.Q)){
          norm.fun <- function(x) {(x - min(x))/(max(x)-min(x)) * (norm[2] - norm[1]) + norm[1]}
          value <- as.data.table(lapply(value, norm.fun)) # by default, removes the names of each row
          
          message("Normalisation between ", norm[1], " and ", norm[2], " complete")
        }
        
        ## With upper Q (upper Q is NOT NULL)
        if(!is.null(upper.Q)){
          for(a in use.cols){
            # a <- use.cols[[8]]
            
            A <- min(value[[a]])
            B <- quantile(value[[a]], probs = upper.Q)
            C <- max(value[[a]])
            
            pt.1 <- (B)-(A) # range from in to upper.Q
            pt.2 <- pt.1 + (C-B) 
            mx <- pt.2/pt.1
            mx <- mx * norm[2]
            
            temp <- value[, a, with = FALSE]
            temp <- as.data.table(temp)
            names(temp) <- a

            # norm.fun <- function(x) {(x - min(x))/(max(x)-min(x)) * (mx - norm[1]) + norm[1]}
            # temp <- as.data.table(norm.fun(temp)) # by default, removes the names of each row
            
            res <- do.normalise(temp, a, new.min = norm[1], new.max = mx)
            res <- res[[paste0(a, "_norm")]]
            
            value[[a]] <- res

            rm(a)
            rm(A)
            rm(B)
            rm(C)
            rm(temp)
            rm(res)
            rm(pt.1)
            rm(pt.2)
          }
          message("Normalisation between ", norm[1], " and ", norm[2], " complete, with upper.Q = ", upper.Q)
        }
      }
      
  ### Wrap up
      
      value <- as.data.table(value)
      names(value) <- paste0(names(value), append.name)
    
      res <- cbind(dat, value)
      
      if(save.info == TRUE){
        dt <- date()
        
        settings <- as.matrix(c(dt,
                                "----------------",
                                paste0("Performed on columns:", use.cols),
                                paste0("Arcsinh co-factor = ", asinh.cf),
                                paste0("Low threshold = ", lower.threshold),
                                paste0("Normalise between = ", norm[1], ' and ', norm[2]),
                                paste0("Upper Q = ", upper.Q),
                                paste0("Appended name = ", append.name)
        ))
        
        write(settings, file = paste0('Multipurpose normaliser (do.mn) settings - ', dt, ".txt"))  
      }
      
      return(res)
}
  
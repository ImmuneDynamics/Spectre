#' do.mn - Multipurpose function for arcsinh transformation, outlier management, and quantile-based re-scaling
#'
#' This function facilitates multiple data transformation functions in a single process, including arcsinh transformation, outlier management, and quantile-based re-scaling.
#'
#' @usage do.mn(dat, use.cols)
#'
#' @param dat NO DEFAULT. A data.table
#' @param use.cols NO DEFAULT. Vector of column names that you wish to transform
#' @param asinh.cf DEFAULT = 15. The co-factor to use for arcsinh transformation
#' @param lower.Q DEFAULT = 0.001. Quantile target for values at the lower threshold. Set to 0 to use the minimum value.
#' @param upper.Q DEFAULT = 0.995. Quantile target for rescaling, such that the data will be scaled between 0 (the minimum value after clipping) and 1 (upper.Q), where values above the upper.Q will increase above 1 proportionatly. Set to 1 to use the maximum value.
#' @param lower.func DEFAULT = 'clip'. Determine what is done to values below the lower.Q threshold. 'clip' means values below the lower.Q quantile threshold will be clipped to that value. NULL will skip.
#' @param append.name DEFAULT = '_transf'. Text to be appended to resultant column name
#' @param save.info DEFAULT = TRUE. Save's the transformation details to a .txt file.
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
#' @import data.table
#'
#' @export

do.mn <- function(dat,
                  use.cols,
                   
                  asinh.cf = 15,
                  
                  lower.Q = 0.001,
                  upper.Q = 0.995,
                  
                  lower.func = 'clip', 
                  
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
      # 
      # lower.Q = 0.005
      # upper.Q = 0.995
      # 
      # lower.func = 'clip'
      # 
      # save.info = TRUE
      # append.name = '_transf'
      # info.dir = getwd()

  ### Checks
  
      if(length(use.cols) == 0){
        return(dat)
      }

      if(lower.Q >= upper.Q){
        stop(paste0("lower.Q is greater than upper.Q"))
      }
      
  ### Extract data
  
      norm = c(0,1) # default re-scaling to 0 and 1
  
      value <- dat[,use.cols,with = FALSE]
  
  ### Check for numeric values
      
      if(isFALSE(all(sapply(value, is.numeric)))){
        message("It appears that one column in your dataset is non numeric")
        print(sapply(value, is.numeric))
        stop("do.asinh stopped")
      }
  
  ### Arcsinh transformation
  
      if(!is.null(asinh.cf)){
        value <- value / asinh.cf
        value <- asinh(value) # value <- log(value + sqrt(value^2 + 1))
        #names(value) <- paste0(names(value), "_asinh")
        
        message("Arcsinh transformation complete")
      }
 
  ### Lower Q clipping

      if(lower.func == 'clip'){
        
        for(a in use.cols){
          # a <- use.cols[[2]]
          
          temp <- value[,a,with = FALSE]
          b <- quantile(temp[[1]], lower.Q)
          
          temp[temp[[a]] < b,] <- b
          value[,a] <- temp
          
          rm(a)
          rm(b)
          rm(temp)
        }
        
        message("Lower quantile clipping is complete")
      }

  ### Rescale between min and upper.Q as 1, with values above upper.Q increasing proportionally above 1
  
      for(a in use.cols){
        # a <- use.cols[[2]]
        
        A <- min(value[[a]])
        B1 <- quantile(value[[a]], probs = lower.Q)
        B2 <- quantile(value[[a]], probs = upper.Q)
        C <- max(value[[a]])
        
        ## Upper Q calculations
        pt.1 <- (B2)-(A)        # range from min to upper.Q
        pt.2 <- pt.1 + (C-B2)   # range from min to upper.Q + range from upper.Q to max (total range)
        mx <- pt.2/pt.1         # ratio of [total range] over [min to upper.Q]
        mx <- mx * norm[2]      # ratio x new target value
        
        ## Select data
        temp <- value[, a, with = FALSE]
        temp <- as.data.table(temp)
        names(temp) <- a
        
        ## Normalise
        norm.fun <- function(x) {(x - min(x))/(max(x)-min(x)) * (mx - norm[1]) + norm[1]}
        res <- as.data.table(lapply(temp, norm.fun))
        
        value[[a]] <- res[[1]]
        
        rm(a)
        rm(A)
        rm(B1)
        rm(B2)
        rm(C)
        rm(temp)
        rm(res)
        rm(pt.1)
        rm(pt.2)
      }
      
      message("Re-scaling between ", norm[1], " and ", norm[2], " complete, with lower.Q = ", lower.Q, " and upper.Q = ", upper.Q)
        
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
                                paste0("Lower Q (squish) = ", lower.Q),
                                paste0("Upper Q (re-scale) = ", upper.Q),
                                paste0("Normalise between = ", norm[1], ' and ', norm[2]),
                                paste0("Appended name = ", append.name)
        ))
        
        write(settings, file = paste0('Multipurpose normaliser (do.mn) settings - ', dt, ".txt"))  
      }
      
      return(res)
}
  
#' do.mn - Multipurpose function for arcsinh transformation, outlier management, and quantile-based re-scaling
#'
#' This function facilitates multiple data transformation functions in a single process, including arcsinh transformation, outlier management, and quantile-based re-scaling.
#'
#' @usage do.mn(dat, use.cols)
#'
#' @param dat NO DEFAULT. A data.table
#' @param use.cols NO DEFAULT. Vector of column names that you wish to transform
#' @param asinh.cf DEFAULT = 15. The co-factor to use for arcsinh transformation
#' @param lower.Q DEFAULT = 0.001. Quantile target for values at the lower threshold
#' @param upper.Q DEFAULT = 0.995. Quantile target for rescaling, such that the data will be scaled between 0 (the minimum value after clipping) and 1 (upper.Q), where values above the upper.Q will increase above 1 proportionatly.
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
      # norm = c(0,1)
      # 
      # lower.Q = 0.005
      # upper.Q = 0.995
      # 
      # lower.func = 'clip'
      # 
      # save.info = TRUE
      # append.name = '_transf'
      # info.dir = getwd()

  ### Some settings and tests

      norm = c(0,1) # default re-scaling to 0 and 1
  
      if(lower.squish >= upper.Q){
        stop(paste0("Qs[1] is greater than Qs[2]"))
      }
      
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
 
  ### Lower Q

      if(!is.null(lower.Q)){
        
        if(!is.null(lower.func)){
        
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
          
        }
      }

  # ### Upper squish
      #     
      #     if(!is.null(upper.squish)){
      #       for(a in use.cols){
      #         
      #         b <- upper.squish
      #         temp <- value[,a,with = FALSE]
      #         temp[temp[[a]] > b,] <- b
      #         value[,a] <- temp
      #         
      #         rm(a)
      #         rm(b)
      #         rm(temp)
      #       }
      #       message("Upper squish clipping is complete")
      #     }
      
  ### Rescale between min and upper Q as 1
      
      if(!is.null(norm)){
        
        ## No upper Q (upper Q is NULL)
            
            if(is.null(upper.Q)){
              norm.fun <- function(x) {(x - min(x))/(max(x)-min(x)) * (norm[2] - norm[1]) + norm[1]}
              value <- as.data.table(lapply(value, norm.fun)) # by default, removes the names of each row
              
              message("Normalisation between ", norm[1], " and ", norm[2], " complete")
            }
        
        ## With lower and upper Q (upper Q is NOT NULL)
        
            if(!is.null(Qs)){
                
                  for(a in use.cols){
                    # a <- use.cols[[2]]
                    
                    A <- min(value[[a]])
                    B1 <- quantile(value[[a]], probs = Qs[1])
                    B2 <- quantile(value[[a]], probs = Qs[2])
                    C <- max(value[[a]])
    
                    ## Upper Q calculations
                    
                        pt.1 <- (B2)-(A)        # range from min to upper.Q
                        pt.2 <- pt.1 + (C-B2)   # range from min to upper.Q + range from upper.Q to max (total range)
                        mx <- pt.2/pt.1         # ratio of [total range] over [min to upper.Q]
                        mx <- mx * norm[2]      # ratio x new target value
    
                    ## Lower Q calculations
                        
                        # mn <- ((A) - (B1)) * (mx / C)
                        # mn
    
                        # pt.1 <- (C)-(B1)        # range from max to lower.Q
                        # pt.2 <- pt.1 + (B1)-(A) # range from max to lower.Q + range from lower.Q to min (total range)
                        # 
                        # x - B1 # range from lowerQ to upper Q (=1)
                        # B1 - A # range from lowerQ to min
                        # 
                        # mn <- (A) / (pt.1)
                        # 
                        # # mn <- pt.2 / (norm[1] - ((B1) - (A)))
                    
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
                    rm(B1)
                    rm(B2)
                    rm(C)
                    rm(temp)
                    rm(res)
                    rm(pt.1)
                    rm(pt.2)
                  }
                
              message("Re-scaling between ", norm[1], " and ", norm[2], " complete, with lower.Q = ", Qs[1], " and upper.Q = ", Qs[2])
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
                                paste0("Lower Q (squish) = ", lower.Q),
                                paste0("Upper Q (re-scale) = ", upper.Q),
                                paste0("Normalise between = ", norm[1], ' and ', norm[2]),
                                paste0("Appended name = ", append.name)
        ))
        
        write(settings, file = paste0('Multipurpose normaliser (do.mn) settings - ', dt, ".txt"))  
      }
      
      return(res)
}
  
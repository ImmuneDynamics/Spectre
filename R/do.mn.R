#' do.mn - Multipurpose function for arcsinh transformation, quantile normalisation, noise reduction, and normalisation
#'
#' This function facilitates multiple data transformation functions in a single process, including arcsinh transformation, noise reduction, and normalisation
#'
#' @usage do.mn(dat, use.cols)
#'
#' @param dat NO DEFAULT. A data.table
#' @param use.cols NO DEFAULT. Vector of column names that you wish to transform
#' @param append.name DEFAULT = '_transf'. Text to be appended to resultant column name
#' @param asinh.cf DEFAULT = 15. The co-factor to use for arcsinh transformation
#' @param divide.by DEFAULT = NULL. Column name to divide data -- where each division is subject to the quantile adjustment, minimum thresholding, and normalisation separately
#' @param Qs DEFAULT = NULL. Vector of minimun and maximum quantiles (out of '1') for thresholding (e.g. c(0.005, 0.995)).
#' @param low.threshold DEFAULT = NULL. Values below this number will be converted to this number.
#' @param norm DEFAULT = c(0,1). New MIN and MAX values for each channel, as c(MIN, MAX). 
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
                  append.name = '_transf',
                   
                  asinh.cf = 15,
                   
                  divide.by = NULL,
                  Qs = NULL, # c(0.005, 0.995)
                  low.threshold = 0,
                  norm = c(0,1),
                   
                  info.dir = getwd()
){
  
  ### Libraries
  
      library(Spectre)
      library(data.table)
  
  ### Test data
      
      # dat <- Spectre::demo.start
      # use.cols <- names(dat)[c(2:10)]
      # append.name = '_transf'
      # 
      # asinh.cf <- 15
      # 
      # divide.by = 'FileName'
      # Qs = c(0.005, 0.995)
      # low.threshold = 0
      # norm = c(0,1)
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
 

  ##############################################################################    
  ################################ No divisions  ###############################     
  ##############################################################################    
          
      if(is.null(divide.by)){ 
        
        ### Quantile threhsolding
        
            if(!is.null(Qs)){
                for(a in use.cols){
                  # a <- use.cols[[1]]
                  
                  lower <- quantile(value[[a]], c(Qs[1], Qs[2]))[1]
                  upper <- quantile(value[[a]], c(Qs[1], Qs[2]))[2]
                  
                  value[[a]][value[[a]] < lower] <- lower
                  value[[a]][value[[a]] > upper] <- upper
                  
                  rm(lower)
                  rm(upper)
                  rm(a)
                }
              message("Quantile adjustment between", Qs[1], " and ", Qs[2], ' complete')
            }
        
        ### Minimum threshold
            
            if(!is.null(low.threshold)){
              for(a in use.cols){
                
                b <- low.threshold
                temp <- value[,a,with = FALSE]
                temp[temp[[a]] < b,] <- b
                value[,a] <- temp
                
                rm(a)
                rm(b)
                rm(temp)
                
              }
              message("Noise reduction complete transformation complete")
            }
    
        ### Normalisation
        
            if(!is.null(norm)){
              norm.fun <- function(x) {(x - min(x))/(max(x)-min(x)) * (norm[2] - norm[1]) + norm[1]}
              value <- as.data.table(lapply(value, norm.fun)) # by default, removes the names of each row
              #names(value) <- paste0(names(value), "_norm")
              
              message("Normalisation between ", norm[1], " and ", norm[2], " complete")
            }
    
      }
          
      
  ##############################################################################    
  ##############################  With divisions  ##############################     
  ##############################################################################    
      
      if(!is.null(divide.by)){ 
        
            all.divs <- value
            rm(value)
            
            all.divs$DIVISION_BARCODE <- c(1:nrow(all.divs))
            all.divs$DIVISION_BY <- dat[[divide.by]]
            divs <- unique(dat[[divide.by]])
            
            res.list <- list()
        
            for(d in divs){
              # d <- divs[[1]]
              
              value <- all.divs[all.divs[['DIVISION_BY']] == d,]
              brcd <- value$DIVISION_BARCODE
              
              value$DIVISION_BY <- NULL
              value$DIVISION_BARCODE <- NULL
            
          ### Quantile threhsolding
              
              if(!is.null(Qs)){
                for(a in use.cols){
                  # a <- use.cols[[1]]
                  
                  lower <- quantile(value[[a]], c(Qs[1], Qs[2]))[1]
                  upper <- quantile(value[[a]], c(Qs[1], Qs[2]))[2]
                  
                  value[[a]][value[[a]] < lower] <- lower
                  value[[a]][value[[a]] > upper] <- upper
                  
                  rm(lower)
                  rm(upper)
                  rm(a)
                }
              }
              
          ### Minimum threshold
              
              if(!is.null(low.threshold)){
                for(a in use.cols){
                  
                  b <- low.threshold
                  temp <- value[,a,with = FALSE]
                  temp[temp[[a]] < b,] <- b
                  value[,a] <- temp
                  
                  rm(a)
                  rm(b)
                  rm(temp)
                  
                }
              }
              
          ### Normalisation
              
              if(!is.null(norm)){
                norm.fun <- function(x) {(x - min(x))/(max(x)-min(x)) * (norm[2] - norm[1]) + norm[1]}
                value <- as.data.table(lapply(value, norm.fun)) # by default, removes the names of each row
                #names(value) <- paste0(names(value), "_norm")
              }
              
              value$POST_BARCODE <- brcd
              res.list[[d]] <- value
              
            }
            
            message("Quantile adjustment between", Qs[1], " and ", Qs[2], ' complete')
            message("Noise reduction complete transformation complete")
            message("Normalisation between ", norm[1], " and ", norm[2], " complete")
            
            res <- rbindlist(res.list)
            res <- setorderv(res, "POST_BARCODE")
            res$POST_BARCODE <- NULL
  
            rm(value)
            value <- res
          }

  ### Wrap up
      
      value <- as.data.table(value)
      names(value) <- paste0(names(value), append.name)
    
      res <- cbind(dat, value)
      
      settings <- as.matrix(c(date(),
                              "----------------",
                              paste0("Arcsinh co-factor: ", asinh.cf),
                              paste0("Quantile min and max: ", Qs[1], ' and ', Qs[2]),
                              paste0("Division for quantile thersholding: ", divide.by),
                              paste0("Low threshold: ", low.threshold),
                              paste0("Normalise between: ", norm[1], ' and ', norm[2])
                              ))

      write(settings, file = 'do.mrn settings.txt')
      
      return(res)
}
  
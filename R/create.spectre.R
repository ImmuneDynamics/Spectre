#' Create Spectre object
#' 
#' @export

create.spectre <- function(data = NULL,
                           analysis = NULL,
                           key = NULL, # or vector
                           name = NULL){
  ###
  
  # data <- list('cyto' = Spectre::demo.start[,-1])
  # key <- NULL
  
  ###
  
      require('Spectre2')
      obj <- new('Spectre')
  
  ###
  
      if(is.null(key)){
        z <- paste0("%0", nchar(nrow(data[[1]])), "d")
        x <- sprintf(z,c(1:nrow(data[[1]])))
        y <- data.table('CellID' = paste0('Cell_', x))
      }
      
      if(!is.null(key)){
        y <- data.table('CellID' = key)
      }
  
  ###
  
      if(!is.null(data)){
        message('Adding data')
        for(i in names(data)){
          obj@data[[i]] <- cbind(y, data[[i]])
          setkey(obj@data[[i]], CellID)
        }
      }
      
      if(!is.null(analysis)){
        message('Adding analysis results')
        for(i in names(analysis)){
          obj@analysis[[i]] <- c(y, analysis[[i]])
          setkey(obj@analysis[[i]], CellID)
        }
      }
  
  ###
      
      obj@key <- 'CellID'
      
  ###
      
      if(is.null(name)){
        obj@name <- ''
      } else {
        obj@name <- name
      }
  
  ###
  
      return(obj)
  
}


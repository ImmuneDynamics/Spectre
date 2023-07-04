#' Create Spectre object
#' 
#' @export

create.spectre <- function(data = NULL,
                           analysis = NULL,
                           key = NULL){
  ###
  
  # data <- list('cyto' = Spectre::demo.start[,-1])
  # key <- NULL
  
  ###
  
  obj <- new('Spectre')
  
  # cellid <- 
  
  z <- paste0("%0", nchar(nrow(data[[1]])), "d")
  x <- sprintf(z,c(1:nrow(data[[1]])))
  y <- data.table('CellID' = paste0('Cell_', x))
  
  ###
  
  if(!is.null(data.frame())){
    for(i in names(data)){
      obj@data[[i]] <- cbind(y, data[[i]])
    }
  }
  
  if(!is.null(analysis)){
    for(i in names(analysis)){
      obj@analysis[[i]] <- c(y, analysis[[i]])
    }
  }
  
  ###
  
  # if(!is.null(key)){
  #   obj@key <- 'CellID'
  # } else {
  #   obj@key <- key
  # }
  
  ###
  
  return(obj)
  
}


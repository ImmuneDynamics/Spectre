#' Spectre object
#'
#' @export

# Spectre <- setClass(Class = "Spectre",
#                     slots = c(
#                       data = "list",
#                       analysis = "list",
#                       other = "list",
#                       key ='vector'
#                     )

Spectre <- setClass(Class = "Spectre",representation(data='list',
                                  analysis="list",
                                  other='list',
                                  key='vector'))


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
      
      if(!is.null(key)){
        obj@key <- 'CellID'
      } else {
        obj@key <- key
      }
  
  ###
      
      return(obj)
      
}



#' @import Spectre
  
setMethod("show",
          "Spectre",
          function(object){

            ###

            require('data.table')

            ###
            s <- object.size(object)/1000/1000/1000
            s <- round(s, 2)

            message(' ')
            message("Spectre OBJECT")
            message(' ')
            message(' - Name  : ', 'Test object')
            # message(' - Date  : ')
            message(' - Size  : ', s, ' Gb')
            message(" - Key   : '", object@key, "'")

            ### @data

            message(' ')
            message("@data")

            for(i in names(object@data)){
              x <- object@data[[i]]
              n <- nrow(x)
              message(' ')
              message('... $', i, ' (', paste0(class(x)[1]), ', nrow: ', nrow(x), ' ncol: ', ncol(x), ')')
              message(' ')
              print(as.data.table(x))
            }

            ### @analysis

            message(' ')
            message("@analysis")

            for(i in names(object@analysis)){
              x <- object@analysis[[i]]
              n <- nrow(x)
              message(' ')
              message('... $', i, ' (', paste0(class(x)[1]), ', nrow: ', nrow(x), ' ncol: ', ncol(x), ')')
              message(' ')
              print(as.data.table(x))
            }

            ### @other
            message(' ')
            message("@other")

            for(i in names(object@analysis)){
              message(' ')
              message('... $', i)
            }

          })
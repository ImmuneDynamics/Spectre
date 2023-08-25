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
                                      key='vector',
                                      name='vector'))

#' @import Spectre
  
setMethod("show",
          "Spectre",
          function(object){

            ###
            
                require('data.table')

            ###
               
                if(length(object@name) > 0){
                  message(' ')
                  message("Spectre object (", object@name, ")")
                } else {
                  message(' ')
                  message("Spectre object")
                }
             
            ### @data
    
                message("   @data")
    
                for(i in names(object@data)){
                  # print(object@data[[i]])
                 
                  message(paste0("      $", i))
                  
                  message('          rows (n= ', formatC(nrow(object@data[[i]]), format = 'd', big.mark = ","), '): ', paste0(dat@data[[i]][[1]][1:5], ', '), ' ...')
                  message('          cols (n= ', ncol(object@data[[i]]), '): ', paste0(names(dat@data[[i]])[2:6], ', '), ' ...')
                }
                
            ### @analysis
    
                message("   @analysis")
    
                for(i in names(object@analysis)){
                  message(paste0("     $", i))
                }

            ### @other
                
                message("   @other")
    
                for(i in names(object@other)){
                  message(paste0("     $", i))
                }

          })
                  
#' spectre object
#'
#' @export

setClass("spectre",representation(data='list',
                                  analysis="list",
                                  other='list',
                                  cellid='vector'))

setMethod("show",
          "spectre",
          function(object){
            
            ###
            
            require(data.table)
            
            ### 
            s <- object.size(object)/1000/1000/1000
            s <- round(s, 2)
            
            message(' ')
            message("SPECTRE OBJECT")
            message(' ')
            message(' - Name  : ', 'Test object')
            # message(' - Date  : ')
            message(' - Size  : ', s, ' Gb')
            message(' - CellID: ', dat@cellid)
            
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
            
            ### 
            # message(' ')
            # message(paste0("cellid =", dat@cellid))
            # message(" ")
            
          })



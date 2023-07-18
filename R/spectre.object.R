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
            # s <- object.size(object)/1000/1000/1000
            # s <- round(s, 2)

            message(' ')
            message("Spectre object")
            message(' ')
            message(" - Name  :", object@name)
            # message(' - Date  : ')
            # message(' - Size  : ', s, ' Gb')
            # message(" - Key   : ", object@key)

            ### @data

            message(' ')
            message("@data")

            for(i in names(object@data)){
              
              # if(ncol(object@data[[i]]) > 11){
              #   b <- ncol(object@data[[i]])
              #   a <- b-4
              #   pre <- object@data[[i]][,1:5]
              #   post <- object@data[[i]][,a:b]
              #   x <- cbind(pre, data.table('...' = '...'), post)
              # } else {
              #   x <- object@data[[i]]
              # }
              # message(' ')
              # message('... $', i, ' (', paste0(class(object@data[[i]])[1]), ', nrow: ', nrow(object@data[[i]]), ' ncol: ', ncol(object@data[[i]]), ')')
              # message(' ')
              # print(as.data.table(x))
              
              cols <- c(1:min(8,ncol(object@data[[i]])))
              
              print(object@data[[i]][,..cols])
              print(object@data[[i]][,1:5])
              print(object@data[[i]][,c('CellID', 'NK11', 'CD3', 'CD45', 'Ly6G')], with = FALSE)
              print(object@data[[i]][c(1:5, 170000000:170000005),])
              print(object@data[[i]])
              
              print(object@data$cyto[,1])
              
            }

            ### @analysis

            message(' ')
            message("@analysis")

            for(i in names(object@analysis)){
              
              if(ncol(object@analysis[[i]]) > 11){
                b <- ncol(object@analysis[[i]])
                a <- b-4
                pre <- object@analysis[[i]][,1:5]
                post <- object@analysis[[i]][,a:b]
                x <- cbind(pre, data.table('...' = '...'), post)
              } else {
                x <- object@analysis[[i]]
              }
              message(' ')
              message('... $', i, ' (', paste0(class(object@analysis[[i]])[1]), ', nrow: ', nrow(object@analysis[[i]]), ' ncol: ', ncol(object@analysis[[i]]), ')')
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
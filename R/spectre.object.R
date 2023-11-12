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

Spectre <- setClass(Class = "Spectre",representation('list',
                                                     name='vector',
                                                     spatial="list",
                                                     other="list"
))

setMethod("show",
          "Spectre",
          function(object){
            
            ### Initial
            
            require('data.table')
            
            if(length(dat@name < 1)){
              message(paste0("Spectre object: '", dat@name, "'"))
            } else {
              message("Spectre object")
            }
            
            ### Datasets
            
            for(i in names(dat)){
              r <- nrow(dat[[i]])
              if(is.null(r)){
                r <- 0
              }
              c <- ncol(dat[[i]])
              if(is.null(c)){
                c <- 0
              }
              message(paste0('  $', i, ' (', r, ' rows, ', c, ' cols)'))
            }
            
            ### Spatial
            
            if(length(dat@spatial)){
              message('  --')
              
              l <- length(names(dat@spatial))

              if(l <= 3){
                for(i in names(dat@spatial)){
                  message(paste0('  @spatial$', i))
                }
              } else {
                for(i in c(1:3)){
                  a <- names(dat@spatial)[i]
                  message(paste0('  @spatial$', a))
                }
                message(paste0('   +', l-3))
              }
            }
            
            ### Other
            
            if(length(dat@other)){
              message('  --')
              
              for(i in names(dat@other)){
                message(paste0('  @other$', i))
              }
            }
            
          })


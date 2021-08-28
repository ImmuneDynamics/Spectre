#' spectre object
#' 
#' @export

### Spatial object

spectre <- setClass(Class = 'spectre', 
                    slots = c(meta = 'data.table',
                              data = 'list', 
                              analysis = 'list',
                              info = 'list'
                    ))

create.spectre <- function(){
  dat <- spectre()
  return(dat)
}

setMethod("show",
          "spectre",
          function(dat) {
            
            ### Setup
            
            require('Spectre')
            require('data.table')
            
            #spacers <- data.table('...  ' = rep('...  ', nrow(dat@meta)))
            #lines <- data.table('       ' = rep('    |', nrow(dat@meta)))
            
            B <- list()
            C <- list()
            
            B.rows <- list()
            B.cols <- list()
            
            C.rows <- list()
            
            ### A
            
            if(length(dat@meta) == 0){
              A <- 'No metadata'
              A.rows <- 0
            } else {
              A <- dat@meta 
              A.rows <- (nrow(dat@meta))
            }
            
            ### B
            
            if(length(dat@data) == 0){
              B <- 'No data'
              B.rows <- 0
            } else {
              
              for(i in names(dat@data)){
                
                if(ncol(dat@data[[i]]) < 6){
                  B[[i]] <- dat@data[[i]]
                  
                } else {
                  B[[i]] <- dat@data[[i]][,c(1:6)]
                  B[[i]][[paste0(' +', ncol(dat@data[[i]]) - 6, ' cols')]] <- '...'
                }
                
                B.rows[[i]] <- nrow(dat@data[[i]])
                B.cols[[i]] <- names(dat@data[[i]])
              }
            }
            
            ### C
            
            if(length(dat@analysis) == 0){
              C <- 'No analysis'
              C.rows <- 0
            } else {
              
              for(i in names(dat@analysis)){
                
                if(ncol(dat@analysis[[i]]) < 6){
                  C[[i]] <- dat@analysis[[i]]
                  
                } else {
                  C[[i]] <- dat@analysis[[i]][,c(1:6)]
                  C[[i]][[paste0(' +', ncol(dat@analysis[[i]]) - 6, ' cols')]] <- '...'
                }
                
                C.rows[[i]] <- nrow(dat@analysis[[i]])
                #C.cols[[i]] <- names(dat@analysis[[i]])
              }
            }
            
            ### Row and column check
            
            unique.nrows <- unique(c(A.rows, unlist(B.rows), unlist(B.rows)))
            
            ### Unique feature names
            
            if(length(dat@data) > 0){
              
              unique.cols <- unique(unlist(B.cols))
              n.features <- length(unique.cols)
              
              if(n.features > 6){
                unique.cols <- c(unique.cols)[c(1:6)] 
                fllw <- paste0(' +', n.features-6, ' ...')
              } else {
                unique.cols <- c(unique.cols)
                fllw <- c(' ')
              }
              
            } else {
              unique.cols <- 'No data columns'
            }
            
            ### Print
            
            message(' ')
            message('SPECTRE OBJECT')
            message(' ')
            
            v <- devtools::package_info('spectre')
            v <- v$ondiskversion
            
            message(c(' - Spectre version: ', v))
            
            if(length(dat@data) > 1){
              message(c(' - Cols: ', n.features, ' detected features (', paste0(unique.cols, " "), fllw, ')'))  
            } else {
              message(c(' - Cols: none'))
            }
            
            if(length(unique.nrows) == 1){
              message(c(' - Rows: ', unique.nrows, ' cells'))
            } else {
              message(c(' - Rows: differing number of rows/cells!'))
            }
            
            message(' ')
            
            message("@meta ============================================================================================")
            message(' ')
            print(A)
            message(' ')
            
            
            message("@data ============================================================================================")
            message(' ')
            
            if(length(dat@data) > 1){
              for(i in names(B)){
                message(c('... $', i, ' (', class(B[[i]])[[1]], ')'))
                message(' ')
                print(B[[i]])
                message(' ')
              }
            } else {
              print(B)
            }
            
            message(' ')
            
            
            message("@analysis ========================================================================================")
            message(' ')
            
            if(length(dat@analysis) > 1){
              for(i in names(C)){
                message(c('... $', i, ' (', class(C[[i]])[[1]], ')'))
                message(' ')
                print(C[[i]])
                message(' ')
              }
            } else {
              print(C)
            }
            
            message(' ')
            
          }
)

#' adjust.files
#' 
#' @param dat List of data.tables, typically one per 'sample'
#' @param colnames.table A data.table of column names. First column contains the current names, second column contains the new names to update to.
#' 
#' @export

adjust.files <- function(dat, colnames.table){
 
  new.names <- colnames.table
  
  for(a in names(dat)){
    # a <- names(dat)[1]
    
    for(i in c(1:nrow(new.names))){
      # i <- 1
      old <- new.names[i,1][[1]]
      new <- new.names[i,2][[1]]
      
      if(length(which(names(dat[[a]]) == old)) == 1){
        names(dat[[a]])[which(names(dat[[a]]) == old)] <- new
      }
    }
  }
  
  return(dat)
}
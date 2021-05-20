#' do.filter.percell
#' 
#' @import data.table
#' 
#' @export

do.filter.percell <- function(spatial.dat,
                              per.cell,
                              to,
                              filter.by,
                              id.col = "ObjectNumber",
                              x.col = "Location_Center_X",
                              y.col = "Location_Center_Y",
                              simplify.cp.colname = TRUE,
                              value.modifier = 65535){
  
  ### Setup
  message("This is a developmental Spectre-spatial function that is still in testing phase with limited documentation. We recommend only using this function if you know what you are doing.")
  
  ### Loop
  
  for(i in names(spatial.dat)){
    # i <- "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac"
    temp <- spatial.dat[[i]]$CPDATA[[per.cell]]
    temp.filtered <- temp[,grepl( filter.by , names(temp)),with = FALSE]
    temp.filtered <- temp.filtered * value.modifier
    temp.filtered <- cbind(ID = temp[[id.col]],x = temp[[x.col]], y = temp[[y.col]], temp.filtered)
    
    if(simplify.cp.colname == TRUE){
      measr <- temp.filtered[,c(1:3),with = FALSE]
      chnls <- temp.filtered[,c(4:length(names(temp.filtered))),with = FALSE]
      
      names(chnls) <- sub('.*\\_c', '', names(chnls))
      temp.names <- names(chnls)
      
      for(b in c(1:length(temp.names))){
        a <- temp.names[[b]]
        if(nchar(a) == 1){
          a <- paste0("0", a)
        }
        temp.names[[b]] <- a
      }
      
      names(chnls) <- temp.names
      neworder <- sort(names(chnls))
      chnls <- setcolorder(chnls, neworder)
      
      temp.filtered <- cbind(measr, chnls)
    }
    
    spatial.dat[[i]]$CPDATA[[to]] <- temp.filtered
  }
  
  ### Return
  return(spatial.dat)
  
}



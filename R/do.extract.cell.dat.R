#' do.extract.cell.dat
#' 
#' @import data.table
#' 
#' @export

do.extract.cell.dat <- function(spatial.dat,
                                target.dat){
  
  ### Setup
  message("This is a developmental Spectre-spatial function that is still in testing phase with limited documentation. We recommend only using this function if you know what you are doing.")
  
  ### Extract
  
  dat.list <- list()
  
  for(i in names(spatial.dat)){
    roi.dat <- spatial.dat[[i]]$CPDATA[[target.dat]]
    nme.vec <- rep(i, nrow(roi.dat))
    
    roi.dat <- cbind("ROI" = nme.vec, roi.dat)
    dat.list[[i]] <- roi.dat
  }
  
  dat.dt <- rbindlist(dat.list, fill = TRUE)
  return(dat.dt)
}





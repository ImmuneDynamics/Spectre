#' do.extract.cell.dat
#'
#' @export

do.extract.cell.dat <- function(dat = spatial.dat,
                                type = "per.cell.means"){

  message("This is a developmental Spectre-spatial function that is still in testing phase with limited documentation. We recommend only using this function if you know what you are doing.")

  extract.list <- list()

  for(i in dat$meta.data$roi.names){
    extract.list[[i]] <- dat$rois[[i]][[type]]
  }

  extracted <- rbindlist(extract.list, fill = TRUE)
  return(extracted)
}



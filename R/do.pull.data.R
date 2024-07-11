#' do.pull.data
#'
#' @param dat NO DEFAULT. Spatial data list.
#' @param target.dat NO DEFAULT. Dataset to pull.
#'
#' @usage do.pull.data(spatial.dat, target.dat)
#' @import data.table
#' 
#' @usage do.pull.data(spatial.dat, target.dat)
#'
#' @export do.pull.data

do.pull.data <- function(spatial.dat,
                         target.dat) {

  ### Setup
  # message("This is a developmental Spectre-spatial function that is still in testing phase with limited documentation. We recommend only using this function if you know what you are doing.")

  require(data.table)

  ### Extract

  dat.list <- list()

  for (i in names(spatial.dat)) {
    roi.dat <- spatial.dat[[i]]@DATA[[target.dat]]
    nme.vec <- rep(i, nrow(roi.dat))

    roi.dat <- cbind("ROI" = nme.vec, roi.dat)
    dat.list[[i]] <- roi.dat
  }

  dat.dt <- rbindlist(dat.list, fill = TRUE)
  return(dat.dt)
}

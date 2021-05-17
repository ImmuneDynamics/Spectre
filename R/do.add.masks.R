#' do.add.masks
#'
#' @param spatial.dat NO DEFAULT. Spatial data list.
#' @param mask.loc NO DEFAULT. Directory of mask files.
#' @param masks NO DEFAULT. Vector of mask file names
#' @param mask.label DEFAULT = "cell_mask"
#' @param mask.ext DEFAULT  = "_mask.tiff"
#' @param correct.extent DEFAULT = TRUE
#' @param flip.y DEFAULT = TRUE
#' @param value.modifier DEFAULT = 65535
#' @param array DEFAULT = FALSE
#'
#' @import data.table
#'
#' @export

do.add.masks <- function(spatial.dat,
                         mask.loc,
                         masks,

                         mask.label = "cell_mask",
                         mask.ext = "_mask.tiff",

                         correct.extent = TRUE,
                         flip.y = TRUE,
                         value.modifier = 65535,
                         array = FALSE
                         ){

  ### Setup
  message("This is a developmental Spectre-spatial function that is still in testing phase with limited documentation. We recommend only using this function if you know what you are doing.")

  spat.names <- names(spatial.dat)
  spat.names

  mask.names <- gsub(mask.ext, "", masks)
  mask.names

  mask.check <- (spat.names == mask.names)
  mask.check

  if(all(mask.check) == FALSE){
    stop('Error -- list of ROIs does not match the list of masks')
  }

  ### Read in mask files

  setwd(mask.loc)

  for(i in spat.names){
    # i <- spat.names[[1]]

    if(array == FALSE){
      mask.img <- readTIFF(paste0(i, mask.ext))
      mask.img <- raster(mask.img)
    }

    if(array == TRUE){
      mask.img <- raster(paste0(i, mask.ext))
    }

    if(correct.extent == TRUE){
      extent(mask.img) <- c(0, dim(mask.img)[2], 0,dim(mask.img)[1]) # Y axis - X axis
    }

    if(flip.y == TRUE){
      mask.img <- flip(mask.img, 'y')
    }

    values(mask.img) <- values(mask.img)*value.modifier
    names(mask.img) <- mask.label

    spatial.dat[[i]]$MASKS[[mask.label]]$maskraster <- mask.img

  }

  message("Returning spatial data object with added masks")
  return(spatial.dat)
}


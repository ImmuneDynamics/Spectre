#' do.add.masks
#'
#' @param spatial.dat NO DEFAULT. Spatial data list.
#' @param mask.loc NO DEFAULT. Directory of mask files.
#' @param masks NO DEFAULT. Vector of mask file names
#' @param mask.label DEFAULT = "cell_mask"
#' @param mask.ext DEFAULT  = "_mask.tif" (will also recognise the '.tiff' file type)
#' @param correct.extent DEFAULT = TRUE
#' @param flip.y DEFAULT = TRUE
#' @param value.modifier DEFAULT = 65535
#' @param array DEFAULT = FALSE
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://github.com/ImmuneDynamics/Spectre}.
#'
#' @import data.table
#'
#' @export

do.add.masks <- function(spatial.dat,
                         mask.loc,
                         masks,

                         mask.label = "cell.mask",
                         mask.ext = "_mask.tif",

                         correct.extent = TRUE,
                         flip.y = TRUE,
                         value.modifier = 65535,
                         array = FALSE
                         ){

  ### Require
  
      require('raster')
      require('rhdf5')
      require('HDF5Array')
  
  ### Setup
      
      message(paste0('Adding masks to ', mask.label))
  
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
        
        message(paste0("  -- processing ", i))
        
        if(grepl('.h5', mask.ext)){
          file.ext <- '.h5'
        } else {
          file.ext <- 'other'
        }
 
        ## If HDF5 files
        
            if(file.ext == '.h5'){
              
              message("     ... reading HDF5 file")
              
              # h5closeAll()
              # h5ls(paste0(i, mask.ext))
              
              mask.img <- h5read(paste0(i, mask.ext), name = 'exported_data')
              mask.img <- array(as.numeric(mask.img), dim(mask.img))
              mask.img <- matrix(mask.img, nrow = dim(mask.img)[3], ncol = dim(mask.img)[2])
              mask.img <- raster(mask.img)
              
            } else {
              
        ## If NOT HDF5 files
              
              if(array == FALSE){
                mask.img <- readTIFF(paste0(i, mask.ext))
                mask.img <- raster(mask.img)
              }
              
              if(array == TRUE){
                
                message("     ... reading array")
                
                mask.img <- raster(paste0(i, mask.ext))
              }

            }

        ## Correct extent
        
            if(correct.extent == TRUE){
              message("     ... correcting extent")
              extent(mask.img) <- c(0, dim(mask.img)[2], 0,dim(mask.img)[1]) # Y axis - X axis
            }
    
        ## Flip Y
        
            if(flip.y == TRUE){
              message("     ... flipping Y")
              mask.img <- flip(mask.img, 'y')
            }
    
        ## Finalise
        
            raster::values(mask.img) <- raster::values(mask.img) * value.modifier
            names(mask.img) <- mask.label
            spatial.dat[[i]]$MASKS[[mask.label]]$maskraster <- mask.img
    
      }
    
      message("Returning spatial data object with added masks")
      return(spatial.dat)
}


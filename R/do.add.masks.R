#' do.add.masks
#'
#' @param dat NO DEFAULT. A list of spatial data objects
#' @param mask.dir NO DEFAULT. Directory of mask files.
#' @param mask.pattern NO DEFAULT. A character pattern that identifies the mask type (e.g. '_mask')
#' @param mask.label NO DEFAULT. What do you want to call the mask
#' @param correct.extent DEFAULT = TRUE
#' @param flip.y DEFAULT = TRUE
#' @param value.modifier DEFAULT = 65535
#' @param HDF5 DEFAULT = FALSE. Can read in HDF5 mask files if desired (advanced use).
#' @param array DEFAULT = FALSE. If the mask TIFF is an array (advanced use).
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://github.com/ImmuneDynamics/Spectre}.
#'
#' @import data.table
#'
#' @export

do.add.masks <- function(dat, # list of spatial objects
                         mask.dir,
                         mask.pattern,
                         mask.label,
                      
                         correct.extent = TRUE,
                         flip.y = TRUE,
                         value.modifier = 65535,
                         HDF5 = FALSE,
                         array = FALSE
                         ){

  ### Require
  
      require('raster')
      require('rhdf5')
      require('HDF5Array')
  
  ###
  
      # dat <- spatial.dat
      # dat
      # mask.dir <- '/Users/thomasa/OneDrive - The University of Sydney (Staff)/Library/Github (public)/Spectre/workflows/Spatial - advanced/data/masks/'
      # mask.pattern <- 'Object Identities'
      # mask.label <- "cell.mask"

      message('Reading in mask files')
  
      if(class(dat) != 'list'){
        stop("Your input data 'dat' needs to be a list of spectre spatial objects")
      }

      if(isTRUE(HDF5)){
        all.files <-list.files(mask.dir, pattern = '.h5')
      } else {
        all.files <- list.files(mask.dir, pattern = '.tif')
      }

      as.matrix(all.files)

  ### Setup
      
      # message(paste0('Adding masks to ', mask.label))
      # 
      # spat.names <- names(spatial.dat)
      # spat.names
      # 
      # mask.names <- gsub(mask.ext, "", masks)
      # mask.names
      # 
      # mask.check <- (spat.names == mask.names)
      # mask.check
      # 
      # if(all(mask.check) == FALSE){
      #   stop('Error -- list of ROIs does not match the list of masks')
      # }

  ### Read in mask files
    
      for(i in names(dat)){
        # i <- 'ROI002' 
        
        fls <- all.files[grepl(i, all.files)]
        fls <- fls[grepl(mask.pattern, fls)]
        fls

        if(substr(mask.dir, nchar(mask.dir), nchar(mask.dir)) != '/'){
          fls <- paste0(mask.dir, '/', fls)
        } else {
          fls <- paste0(mask.dir, fls)
        }
        
        ## Read in mask files
        
        message("  -- reading mask '", mask.pattern, "' for ROI", i)

        if(isTRUE(HDF5)){
          message("     ...reading HDF5 file")
          
          mask.img <- h5read(fls, name = 'exported_data')
          mask.img <- array(as.numeric(mask.img), dim(mask.img))
          mask.img <- matrix(mask.img, nrow = dim(mask.img)[3], ncol = dim(mask.img)[2])
          mask.img <- raster(mask.img)
          
        } else {
          
          if(array == FALSE){
            message("     ...reading image file")
            mask.img <- readTIFF(fls)
            mask.img <- raster(mask.img)
          }
          
          if(array == TRUE){
            message("     ... reading array")
            mask.img <- raster(fls)
          }
          
        }
        
        ## Correct extent
        
        if(correct.extent == TRUE){
          message("     ...correcting extent")
          extent(mask.img) <- c(0, dim(mask.img)[2], 0,dim(mask.img)[1]) # Y axis - X axis
        }
        
        ## Flip Y
        
        if(flip.y == TRUE){
          message("     ...flipping Y")
          mask.img <- flip(mask.img, 'y')
        }
        
        ## Finalise
        
        raster::values(mask.img) <- raster::values(mask.img) * value.modifier
        names(mask.img) <- mask.label
        
        spatial.dat[[i]]@MASKS[[mask.label]]$maskraster <- mask.img
      }
      
  ### Return
      
      return(spatial.dat)
      
}


#' Read TIFF files into R and and create spatial data object
#'
#' @usage read.spatial.files()
#'
#' @param dir Directory location of where the ROI files are located
#' @param lois Default = NULL. You can provide a list of only some of the ROI folders to read in.
#' @param correct.extent Correct extent so that the minimum is 0,0. Default TRUE
#' @param flip.y Flip the data arrangement for the y-axis. Default TRUE
#' @param value.modifier Data modifier based on image processing. Default 65535
#'
#' @return Returns a list of spatial data objects.
#'
#' @examples
#'
#' @import data.table
#'
#' @export

read.spatial.files <- function(dir,
                               rois = NULL,
                               correct.extent = TRUE,
                               flip.y = TRUE,
                               value.modifier = 65535
                               ){

  ### Packages
  
      require('data.table')
      require('raster')
      require('raster')
      require('tiff')

  ### Setup
  
      if(is.null(rois)){
        rois <- list.files(dir)
      }

      spatial.dat <- spatial()
      
      message('Reading TIFFs from:', dir)
      
      if(correct.extent == TRUE){
        message("  ...with extent correction")
      }
      
      if(flip.y == TRUE){
        message("  ...with y-axis orientation flipping")
      }
      
  ### Read in files
      
      roi.list <- list()
      
      for(a in rois){
        # a <- rois[[1]]
        
        message("Reading ROI: ", a)
        
        if(substr(dir, nchar(dir), nchar(dir)) != '/'){
          fls <- list.files(paste0(dir, '/', a))
        } else {
          fls <- list.files(paste0(dir, a))
        }

        tiff.list <- list()
        
        for(i in fls){
          # i <- fls[[1]]
          
          ext <- tools::file_ext(i)
          ext <- paste0('.', ext)
          
          tiff.list[[i]] <- tiff::readTIFF(paste0(dir, '/', a, '/', i))
          
          if('matrix' %in% class(tiff.list[[i]])){
            message("  -- reading single band in TIFF:", i)
            tiff.list[[i]] <- raster(tiff.list[[i]])
          }
          if(! is.matrix(tiff.list[[i]])){
          # if(class(tiff.list[[i]]) != 'matrix'){
            if('array' %in% class(tiff.list[[i]])){

              message("  -- merging multiple bands in TIFF:", i)

              tiff.list[[i]] <- NULL
              tiff.list[[i]] <- raster(paste0(dir, '/', a, '/', i))
              nbands <- tiff.list[[i]]@file@nbands
              # nbands

              band.list <- list()

              for(u in c(1:nbands)){
                # band.list[[u]] <- raster::values(raster(tiff.list[[i]], band = u))
                band.list[[u]] <- raster::values(raster(paste0(dir, '/', a, '/', i), band = u))
              }

              band.res <- band.list[[1]]

              for(u in c(2:nbands)){
                band.res <- band.res + band.list[[u]]
              }

              # tiff.list[[i]] <- NULL
              # tiff.list[[i]] <- raster(paste0(dir, '/', a, '/', i))
              raster::values(tiff.list[[i]]) <- band.res

            }
          }
          
          if(correct.extent == TRUE){
            extent(tiff.list[[i]]) <- c(0, dim(tiff.list[[i]])[2], 0,dim(tiff.list[[i]])[1]) # Y axis - X axis
          }
          
          if(flip.y == TRUE){
            tiff.list[[i]] <- flip(tiff.list[[i]], 'y')
          }
          
          raster::values(tiff.list[[i]]) <- raster::values(tiff.list[[i]]) * value.modifier
          
          for(n in c(1:length(names(tiff.list)))){
            names(tiff.list)[n] <- gsub(ext, "", names(tiff.list)[n])
          }
          
        }
        
        raster.stack <- stack(tiff.list)
        spatial.dat@RASTERS <- raster.stack
        roi.list[[a]] <- spatial.dat
        
      }
      
      return(roi.list)

}

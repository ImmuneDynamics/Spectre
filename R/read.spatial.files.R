#' Read TIFF files into R and and create spatial data object
#'
#' @usage read.spatial.files()
#'
#' @param rois vector of ROI names (directory names)
#' @param roi.loc working directory for ROIs
#' @param multi.tiff Default FALSE. Currently we don't support reading in multipage TIFF files
#' @param correct.extent Correct extent so that the minimum is 0,0. Default TRUE
#' @param flip.y Flip the data arrangement for the y-axis. Default TRUE
#' @param value.modifier Data modifier based on image processing. Default 65535
#' @param ext Default = ".tif" (which will recognise '.tiff' files as well)
#'
#' @return Returns a spatial data object.
#'
#' @examples
#'
#' @import data.table
#'
#' @export

read.spatial.files <- function(roi.dir,
                               files = NULL,
                               correct.extent = TRUE,
                               flip.y = TRUE,
                               value.modifier = 65535
                               #ext = ".tif"
                               ){

  ### Packages
  
      require('raster')
      require('tiff')

  ### Setup
  
      # roi.dir <- '~/OneDrive - The University of Sydney (Staff)/Library/Github (public)/Spectre/workflows/Spatial - advanced/data/ROIs/ROI002'
      # roi.dir

      if(is.null(files)){
        files <- list.files(roi.dir)
      }

      spatial.dat <- spatial()
      
      message('Reading TIFFs from:', roi.dir)
      
      if(correct.extent == TRUE){
        message("  ...with extent correction")
      }
      
      if(flip.y == TRUE){
        message("  ...with y-axis orientation flipping")
      }
      
  ### Read in files
      
      file.list <- list()
      
      for(a in files){
        # a <- files[[1]]
        
        ext <- tools::file_ext(a)
        ext <- paste0('.', ext)
        
        if(substr(roi.dir, nchar(roi.dir), nchar(roi.dir)) != '/'){
          file.list[[a]] <- tiff::readTIFF(paste0(roi.dir, '/', a))
        } else {
          file.list[[a]] <- tiff::readTIFF(paste0(roi.dir, a))
        }
        
        
        if('matrix' %in% class(file.list[[a]])){
          
          message("  -- reading single band in TIFF:", a)
          file.list[[a]] <- raster(file.list[[a]])
        }
        
        
        if(class(file.list[[a]]) != 'matrix'){
          if('array' %in% class(file.list[[a]])){
            
            message("  -- merging multiple bands in TIFF:", a)
            
            file.list[[a]] <- raster(a)
            
            nbands <- file.list[[a]]@file@nbands
            # nbands
            
            band.list <- list()
            
            for(u in c(1:nbands)){
              band.list[[u]] <- values(raster(a, band = u))
            }
            
            band.res <- band.list[[1]]
            
            for(u in c(2:nbands)){
              band.res <- band.res + band.list[[u]]
            }
            
            values(file.list[[a]]) <- band.res
            
          }
        }
        
        if(correct.extent == TRUE){
          extent(file.list[[a]]) <- c(0, dim(file.list[[a]])[2], 0,dim(file.list[[a]])[1]) # Y axis - X axis
        }
        
        if(flip.y == TRUE){
          file.list[[a]] <- flip(file.list[[a]], 'y')
        }
        
        raster::values(file.list[[a]]) <- raster::values(file.list[[a]]) * value.modifier
        
        for(n in c(1:length(names(file.list)))){
          names(file.list)[n] <- gsub(ext, "", names(file.list)[n])
        }
      }
      
  ### Condense    
  
      raster.stack <- stack(file.list)
      spatial.dat@RASTERS <- raster.stack

      return(spatial.dat)

}

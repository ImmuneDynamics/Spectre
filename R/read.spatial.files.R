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
#' @param ext Default = ".tiff". Could be '.tif' depending on the files.
#'
#' @return Returns a spatial data object.
#'
#' @examples
#'
#' @import data.table
#'
#' @export

read.spatial.files <- function(rois,
                               roi.loc = getwd(),
                               multi.tiff = FALSE,
                               correct.extent = TRUE,
                               flip.y = TRUE,
                               value.modifier = 65535,
                               ext = ".tiff"){

  ### Packages
  
      require('raster')
  
  ### Checks

      if(multi.tiff == TRUE){
        if(length(grep(".tif", rois)) != length(rois)){
          stop("It appears that your list of ROIs are not TIFF stack files, and might be directories full of single TIFFs (i.e. one TIFF per channel. If this is correct, please use 'multi.tiff = FALSE'")
        }
      }

  ### Setup
  
      # message("This is a developmental Spectre-spatial function that is still in testing phase with limited documentation. We recommend only using this function if you know what you are doing.")
      ROI.list <- list()
      spatial.dat <- list()
      setwd(roi.loc)

  ### Loop for ROIs -- one TIFF per ROI (i.e. tiff stack)

      if(multi.tiff == TRUE){
        message("Multi.tiff is not currently supported")
        #
        #   setwd(roi.loc)
        #   ROI.list <- list()
        #
        #   for(i in rois){
        #     # i <- rois[[1]]
        #     message(paste0("Reading TIFF file for ", i))
        #
        #     temp <- readTIFF(i, all = TRUE)
        #     str(temp)
        #
        #     temp <- raster(temp[[1]])
        #
        #     raster.stack <- stack(active.roi)
        #
        #
        #
        #
        #
        #     raster.stack <- stack(active.roi)
        #     ROI.list[[i]]$rasters <- raster.stack
        #   }
      }

  ### Loop for ROIs -- one FOLDER per ROI

      if(multi.tiff == FALSE){
        for(i in rois){
          # i <- rois[[1]]
          message(paste0("Reading TIFF files for ", i))

          setwd(roi.loc)
          setwd(i)
          tiffs <- list.files(pattern = ext)

          ## TIFF loop
          active.roi <- list()

          for(a in tiffs){
            # a <- tiffs[[7]]
            
            active.roi[[a]] <- readTIFF(a)
            
            if('matrix' %in% class(active.roi[[a]])){
              
              message("  -- ", a, " - single band in TIFF")
              active.roi[[a]] <- raster(active.roi[[a]])
              
            }
            
            
            if(class(active.roi[[a]]) != 'matrix'){
              if('array' %in% class(active.roi[[a]])){
                
                message("  -- ", a, " - merging multiple bands in TIFF ")
                
                active.roi[[a]] <- raster(a)
                
                nbands <- active.roi[[a]]@file@nbands
                # nbands
                
                band.list <- list()
  
                for(u in c(1:nbands)){
                  band.list[[u]] <- values(raster(a, band = u))
                }
                
                band.res <- band.list[[1]]
       
                for(u in c(2:nbands)){
                  band.res <- band.res + band.list[[u]]
                }
                
                values(active.roi[[a]]) <- band.res
              
              }
            }
            
            if(correct.extent == TRUE){
              message("    ...correcting extent")
              extent(active.roi[[a]]) <- c(0, dim(active.roi[[a]])[2], 0,dim(active.roi[[a]])[1]) # Y axis - X axis
            }

            if(flip.y == TRUE){
              message("    ...flipping y-axis orientation")
              active.roi[[a]] <- flip(active.roi[[a]], 'y')
            }

            raster::values(active.roi[[a]]) <- raster::values(active.roi[[a]]) * value.modifier

            for(n in c(1:length(names(active.roi)))){
              names(active.roi)[n] <- gsub(ext, "", names(active.roi)[n])
            }
          }

          raster.stack <- stack(active.roi)
          ROI.list[[i]]$RASTERS <- raster.stack
        }
      }


  ### Return
      message("Spatial data object construction complete")
      return(ROI.list)

}

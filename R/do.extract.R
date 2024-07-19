#' do.extract
#'
#' @param dat NO DEFAULT. List of spatial data objects
#' @param mask NO DEFAULT. Name of mask to use for cell data extraction
#' @param name DEFAULT = 'CellData'. Name of the extracted cell dataset.
#' @param fun DEFAULT = 'mean'. Summarisation function.
#'
#' @import data.table
#' 
#' @usage do.extract(dat, mask, name = "CellData", fun = "mean")
#'
#' @export do.extract

do.extract <- function(dat, # spatial.data object
                       mask, # name of the mask being summarised
                       name = "CellData",
                       fun = "mean" # type of marker summarisation (mean, median etc)
                       ){

  #message("This is a developmental Spectre-spatial function that is still in testing phase with limited documentation. We recommend only using this function if you know what you are doing.")

      require('sp')
      require('sf')
      require('exactextractr')
      require('data.table')

  ### Demo data

      # spatial.dat <- demo.spatial
      # str(spatial.dat, 3)
      # 
      # for(i in names(spatial.dat)){
      #   spatial.dat[[i]]$DATA <- NULL
      # }
      # 
      # str(spatial.dat, 3)
      # 
      # dat <- spatial.dat
      # mask <- "cell.mask"
      # name = "CellData"
      # fun = "mean"
      # 
      # str(dat, 4)

  ### Loop for each ROI

      rois <- names(dat)

      for(roi in rois){
        # roi <- rois[[1]]
        
        message(paste0("Processing ", paste0(which(roi == rois), ' of ', length(rois)), ':  ', roi))
 
        roi.stack <- dat[[roi]]@RASTERS
        roi.poly <- dat[[roi]]@MASKS[[mask]]$polygons

        raster.names <- names(roi.stack)
        ply.df <- as.data.frame(roi.poly)
        ply.df

        ply.centroids <- sp::coordinates(roi.poly)
        
        ply.centroids.df <- as.data.frame(ply.centroids)
        ply.centroids.df # mask number, with X and Y coordinates

        # ply.centroids.df <- cbind(ply.centroids.df, as.data.frame(area(roi.poly)))
        
        res <- data.table()
        
        for(i in c(1:length(roi.poly@polygons))){
          x <- data.table('ID' = i, 'Area' = roi.poly@polygons[[i]]@area)
          res <- rbindlist(list(res, x))
          rm(x)
        }
        
        ply.centroids.df <- cbind(ply.centroids.df, res[,2])
        
        names(ply.centroids.df) <- c("x", "y", "Area")

        ## RASTERS

            message(" -- running rasters")
        
            for(i in raster.names){
              # i <- raster.names[[1]]
              message(paste0("     ... ", i))
              temp.dat <- roi.stack[[i]]

              ## Slower method

                  # extracted.dat <- raster::extract(x = temp.dat, y = roi.poly, df = TRUE) # this is the time consuming step
                  # extracted.dat.res <- aggregate(. ~ID, data = extracted.dat, FUN = fun)
                  # # #colnames(extracted.dat.res)[2] <- i # should we be removing .tiff here? If we do should be the same in the other read.spatial function, to ensure matching consistency
                  #
                  # ply.centroids.df <- cbind(ply.centroids.df, extracted.dat.res[2]) ## doing this would remove the necessity to calculate centroids within the 'make.spatial.plot' function

              ## FAST method
              
                  ## Faster options

                  # vx <- velox(temp.dat)
                  # res <- vx$extract(sp=roi.poly, fun=mean) # 3293 polygons #3294?
                  
                  res <- exactextractr::exact_extract(temp.dat, st_as_sf(roi.poly), fun=fun)
                  res <- as.data.table(res)
                  
                  names(res) <- i
                  ply.centroids.df <- cbind(ply.centroids.df, res) ## doing this would remove the necessity to calculate centroids within the 'make.spatial.plot' function
            }

            ID <- c(1:nrow(ply.centroids.df))
            roi.dat <- cbind(ID, ply.centroids.df)
            roi.dat <- as.data.table(roi.dat)

        ## OTHER MASK POLYGONS

            message(" -- running mask polygons")
            
            other.polys <- names(dat[[roi]]@MASKS)
            other.polys <- other.polys[!other.polys %in% mask]

            cols <- c("x", "y", "ID")

            roi.dat.xyid <- roi.dat[,..cols]
            names(roi.dat.xyid) <- c('Longitude', 'Latitude', 'Names')
            roi.dat.xyid

            Longitude <- roi.dat.xyid$Longitude
            Latitude <- roi.dat.xyid$Latitude
            coordinates(roi.dat.xyid) <- ~ Longitude + Latitude

            if(length(other.polys) != 0){
              for(i in c(1:(length(other.polys)))){
                # i <- 1
                
                ply.name <- other.polys[[i]]

                message(paste0("     ... occurance in ", ply.name))

                ply <- dat[[roi]]@MASKS[[ply.name]]$polygons

                proj4string(roi.dat.xyid) <- proj4string(ply)

                over.res <- over(roi.dat.xyid, ply)
                over.res <- as.data.table(over.res)
                roi.dat <- cbind(roi.dat, over.res)
              }
            }

            dat[[roi]]@DATA[[name]] <- roi.dat
      }

  ### Return new dat object
      return(dat)
}

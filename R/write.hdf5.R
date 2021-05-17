#' write.hdf5
#'
#' @param dat NO DEFAULT. Spatial data list
#' @param channels DEFAULT = NULL. Channels to include. If NULL, all channels included.
#' @param q.max DEFAULT = 0.99
#' @param flip.y DEFAULT = TRUE
#' @param random.crop.x DEFAULT = NULL
#' @param random.crop.y DEFAULT = NULL
#' @param random.crop.x.seed DEFAULT = 42
#' @param random.crop.y.seed DEFAULT = 21
#' @param print.spatial.plot DEFAULT = TRUE
#' @param chunk.size DEFAULT = NULL
#' @param compression DEFAULT = 0
#' @param plots DEFAULT = TRUE
#' 
#' @import data.table
#'
#' @export

# run.spatial.analysis

write.hdf5 <- function(dat, # SpectreMAP object
                       channels = NULL,
                       q.max = 0.99,
                       flip.y = TRUE,
                       random.crop.x = NULL,
                       random.crop.y = NULL,
                       random.crop.x.seed = 42,
                       random.crop.y.seed = 21,
                       print.spatial.plot = TRUE,
                       chunk.size = NULL,
                       compression = 0,
                       plots = TRUE
                       ){

  ### Packages
  
      require('Spectre')
      require('data.table')
      require('raster')
      require('rhdf5')
      require('HDF5Array')
  
  ### Testing
      
      # require('SpectreMAP')
      # dat <- spatial.dat
      # channels = names(spatial.dat[[1]]$RASTERS)[c(13:25)]
      # q.max = 0.99
      # flip.y = TRUE
      # random.crop.x = 300
      # random.crop.y = 300
      # random.crop.x.seed = 42
      # random.crop.y.seed = 21
      # print.spatial.plot = TRUE
      # compression = 0
  
  ### Setup

      message('Writing HDF5 files to disk...')

  ### Loop
  
      for(i in names(dat)){
        # i <- names(dat)[1]
        
        ## Flip Y
        
            if(isTRUE(flip.y)){
              dat[[i]]$RASTERS <- raster::flip(dat[[i]]$RASTERS, direction='y')
            }

        ## Filtering to clip very bright pixels
        
            for(a in channels){
              # a <- for.ilastik[1]
              
              max(dat[[i]]$RASTERS[[a]]@data@values)
              
              q <- quantile(dat[[i]]$RASTERS[[a]]@data@values, q.max)
              
              raster::values(dat[[i]]$RASTERS[[a]])[dat[[i]]$RASTERS[[a]]@data@values > q] <- q

              max(dat[[i]]$RASTERS[[a]]@data@values)
              
              if(max(raster::values(dat[[i]]$RASTERS[[a]])) > q){
                stop('Upper threshold clipping did not work')
              }
              
              if(max(dat[[i]]$RASTERS[[a]]@data@values) > q){
                stop('Upper threshold clipping did not work')
              }
              
            }
        
        ## Cropping
            
            if(!is.null(random.crop.x)){
              if(!is.null(random.crop.y)){
                
                ## Dimensions
                    
                    xmin <- extent(dat[[i]]$RASTERS)[1]
                    xmax <- extent(dat[[i]]$RASTERS)[2]
                    ymin <- extent(dat[[i]]$RASTERS)[3]
                    ymax <- extent(dat[[i]]$RASTERS)[4]
                
                ## New X
                
                    x.length <- xmax - xmin
                    x.mid <- xmax - random.crop.x
                    
                    if(x.mid < xmin){
                      x.mid <- xmin    
                    }
                    
                    set.seed(random.crop.x.seed)
                    new.xmin <- runif(1, xmin, x.mid)
                    new.xmin <- ceiling(new.xmin)
                    
                    new.xmax <- new.xmin + random.crop.x
                    
                    if(new.xmax > xmax){
                      new.xmax <- xmax 
                    }
                    
                ## New Y
                
                    y.length <- ymax - ymin
                    y.mid <- ymax - random.crop.y
                    
                    if(y.mid < ymin){
                      y.mid <- ymin    
                    }
                    
                    set.seed(random.crop.y.seed)
                    new.ymin <- runif(1, ymin, y.mid)
                    new.ymin <- ceiling(new.ymin)
                    
                    new.ymax <- new.ymin + random.crop.y 
                    
                    if(new.ymax > ymax){
                      new.ymax <- ymax 
                    }
                
                ## Crop

                    e <- extent(new.xmin, new.xmax, new.ymin, new.ymax)
                    dat[[i]]$RASTERS <- crop(dat[[i]]$RASTERS, e)
              }
            }

        ## Adjust
            
            rws <- dat[[i]]$RASTERS@nrows
            cls <- dat[[i]]$RASTERS@ncols
            
            temp <- raster::values(dat[[i]]$RASTERS)
            temp <- as.data.table(temp)
            temp <- temp[,..channels]
            temp <- as.vector(as.matrix(temp))
            
        ## Array
            
            AR <- array(data = temp, dim = c(rws, cls, length(for.ilastik), 1))
            AR <- aperm(AR, c(3,1,2,4))
        
            dim(AR)[1] # layers (channels)
            dim(AR)[2] # rows (y-axis)
            dim(AR)[3] # columns (x-axis)
            dim(AR)[4] # 1 (z-axis, required for easy reading by Ilastik)
            
        ## Write HDF5
            
            file.remove(paste0(i, ".h5"))
            h5createFile(paste0(i, ".h5"))
            h5createDataset(paste0(i, ".h5"), 
                            "stacked_channels", 
                            dim(AR),
                            storage.mode = "integer", 
                            chunk = c(1,dim(AR)[2],dim(AR)[3],1),
                            level = compression
                            )
            h5write(AR, file=paste0(i, ".h5"),
                    name="stacked_channels")
            h5ls(paste0(i, ".h5"))
            
        ## Message
            
            message('  -- HDF5 file for ', i, ' complete')

        ## Plots
            
            if(isTRUE(plots)){
              
              dir.create('HDF5 plots')
              dir.create(paste0('HDF5 plots/', i))
              
              for(a in channels){
                make.spatial.plot(dat,
                                  image.roi = i,
                                  image.channel = a,
                                  image.min.threshold = 0,
                                  image.max.threshold = 1,
                                  image.y.flip = FALSE, 
                                  path = paste0('HDF5 plots/', i)
                )
              }
              
            }

      }
 
}
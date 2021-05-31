#' write.hdf5
#'
#' @param dat NO DEFAULT. Spatial data list
#' @param channels DEFAULT = NULL. Channels to include. If NULL, all channels included.
#' @param merge.channels DEFAULT = NULL. Channels to use to create an additional 'merged' channel.
#' @param q.max DEFAULT = 0.99
#' @param flip.y DEFAULT = TRUE
#' @param random.crop.x DEFAULT = NULL
#' @param random.crop.y DEFAULT = NULL
#' @param random.crop.x.seed DEFAULT = 42
#' @param random.crop.y.seed DEFAULT = 21
#' @param print.spatial.plot DEFAULT = TRUE
#' @param chunk.size DEFAULT = NULL
#' @param compression DEFAULT = 0
#' @param plots DEFAULT = FALSE
#' 
#' @import data.table
#'
#' @export

# run.spatial.analysis

write.hdf5 <- function(dat, # SpectreMAP object
                       channels = NULL,
                       merge.channels = NULL,
                       q.max = 0.99,
                       flip.y = TRUE,
                       random.crop.x = NULL,
                       random.crop.y = NULL,
                       random.crop.x.seed = 42,
                       random.crop.y.seed = 21,
                       print.spatial.plot = TRUE,
                       chunk.size = NULL,
                       compression = 0,
                       plots = FALSE
                       ){

  ### Packages
  
      require('Spectre')
      require('data.table')
      require('raster')
      require('rhdf5')
      require('HDF5Array')
  
  ### Testing
  
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
      message('...')

  ### Loop
      
      for(i in names(dat)){
        # i <- names(dat)[1]
        
        gc()
        
        roi.dat <- list()
        roi.dat[[i]]<- dat[[i]]
        
        message('Starting ', i)
        
        ## Flip Y
        
            if(isTRUE(flip.y)){
              message("  -- flipping y-axis")
              roi.dat[[i]]$RASTERS <- raster::flip(roi.dat[[i]]$RASTERS, direction='y')
            }

        ## Filtering to clip very bright pixels
        
            message("  -- clipping bright pixels")
            
            for(a in channels){
              # a <- for.ilastik[1]

              max(roi.dat[[i]]$RASTERS[[a]]@data@values)
              
              q <- quantile(roi.dat[[i]]$RASTERS[[a]]@data@values, q.max)
              
              raster::values(roi.dat[[i]]$RASTERS[[a]])[roi.dat[[i]]$RASTERS[[a]]@data@values > q] <- q

              max(roi.dat[[i]]$RASTERS[[a]]@data@values)
              
              if(max(raster::values(roi.dat[[i]]$RASTERS[[a]])) > q){
                stop('Upper threshold clipping did not work')
              }
              
              if(max(roi.dat[[i]]$RASTERS[[a]]@data@values) > q){
                stop('Upper threshold clipping did not work')
              }
            }
        
        ## Cropping
            
            if(!is.null(random.crop.x)){
              if(!is.null(random.crop.y)){
                
                ## Dimensions
                
                    message("  -- setup cropping")
                    
                    xmin <- extent(roi.dat[[i]]$RASTERS)[1]
                    xmax <- extent(roi.dat[[i]]$RASTERS)[2]
                    ymin <- extent(roi.dat[[i]]$RASTERS)[3]
                    ymax <- extent(roi.dat[[i]]$RASTERS)[4]
                
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

                    message("  -- calculate new extent")
                    e <- extent(new.xmin, new.xmax, new.ymin, new.ymax)
                    
                    gc()
                    message("  -- cropping ROI")
                    roi.dat[[i]]$RASTERS <- crop(roi.dat[[i]]$RASTERS, e)
                    
                ## Clean up
                    
                    rm(xmin)
                    rm(xmax)
                    rm(ymin)
                    rm(ymax)
                    rm(e)
                    
                    gc()
              }
            }

        ## Adjust
            
            message("  -- adjusting data")
        
            rws <- roi.dat[[i]]$RASTERS@nrows # rows = y-axis
            cls <- roi.dat[[i]]$RASTERS@ncols # cols = x-axis
            
            temp <- raster::values(roi.dat[[i]]$RASTERS)
            temp <- as.data.table(temp)
            temp <- temp[,..channels]
            
            if(!is.null(merge.channels)){
              
              # sapply(temp[,..merge.channels], max, na.rm = TRUE)
              # colSums(temp[,..merge.channels])
              # rw.tbl <- temp[,..merge.channels]
              
              rw.tbl <- do.rescale(temp, merge.channels, new.min = 0, new.max = 10)
              
              # rw.tbl <- do.rescale(temp, merge.channels, new.min = 0, new.max = max(temp[,..merge.channels]))
              rw.tbl <- rw.tbl[,names(rw.tbl)[grepl('_rescaled', names(rw.tbl))], with = FALSE]

              for(u in names(rw.tbl)){
                # u <- names(rw.tbl)[[5]]

                # rw.tbl[[u]][rw.tbl[[u]] > quantile(rw.tbl[[u]], 0.95)] <- quantile(rw.tbl[[u]], 0.95)
                over.one <- rw.tbl[[u]][rw.tbl[[u]] > 0]
                rw.tbl[[u]][rw.tbl[[u]] <= quantile(over.one, 0.05)] <- 0
              }

              # sapply(rw.tbl, max, na.rm = TRUE)
              # rw.sms <- rowMeans(rw.tbl)

              rw.sms <- rowSums(rw.tbl)
              # rw.sms[rw.sms > length(merge.channels)] <- length(merge.channels)
              # rw.sms[rw.sms > quantile(rw.sms, 0.95)] <- quantile(rw.sms, 0.95)

              
              
              # rw.tbl <- temp[,..merge.channels]
              # rw.tbl <- do.zscore(rw.tbl, merge.channels, replace = TRUE)
              # 
              # for(u in names(rw.tbl)){
              #   # u <- names(rw.tbl)[[5]]
              #   rw.tbl[[u]][rw.tbl[[u]] > quantile(rw.tbl[[u]], 0.95)] <- quantile(rw.tbl[[u]], 0.95)
              #   rw.tbl[[u]][rw.tbl[[u]] <= quantile(rw.tbl[[u]], 0.05)] <- quantile(rw.tbl[[u]], 0.05)
              # }
              
              # rw.sms <- rowSums(rw.tbl)
              rw.sms[rw.sms > quantile(rw.sms, 0.95)] <- quantile(rw.sms, 0.95)
              rw.sms[rw.sms <= quantile(rw.sms, 0.05)] <- 0
              
              temp$Merged <- rw.sms
              temp <- temp[,c('Merged', channels), with = FALSE]
            }
            
            lgth <- length(names(temp))
            temp <- as.vector(as.matrix(temp))
            
        ## Array
            
            message("  -- creating array")
            
            AR <- array(data = temp, dim = c(cls, rws, lgth, 1))
            AR <- aperm(AR, c(3,1,2,4))
        
            rm(temp)
            gc()
            
            dim(AR)[1] # layers (channels)
            dim(AR)[2] # rows (y-axis)
            dim(AR)[3] # columns (x-axis)
            dim(AR)[4] # 1 (z-axis, required for easy reading by Ilastik)
            
        ## Write HDF5
            
            # h5closeAll()
            
            message("  -- writing HDF5 file")
            
            if(!is.null(random.crop.x)){
              if(!is.null(random.crop.y)){
                
                if(!paste0(i, "_crop.h5") %in% list.files(getwd())){
                  file.remove(paste0(i, "_crop.h5"))
                }

                h5createFile(paste0(i, "_crop.h5"))
                h5createDataset(paste0(i, "_crop.h5"), 
                                "stacked_channels", 
                                dim(AR),
                                storage.mode = "integer", 
                                chunk = c(1, dim(AR)[2],dim(AR)[3], 1),
                                level = compression
                                )
                h5write(AR, file=paste0(i, "_crop.h5"),
                        name="stacked_channels")
                h5ls(paste0(i, "_crop.h5"))
                
                  
              }
              
            } else {
              
              if(!paste0(i, ".h5") %in% list.files(getwd())){
                file.remove(paste0(i, ".h5"))
              }
              
              h5createFile(paste0(i, ".h5"))
              h5createDataset(paste0(i, ".h5"), 
                              "stacked_channels", 
                              dim(AR),
                              maxdims = dim(AR),
                              storage.mode = "integer", 
                              chunk = c(1, dim(AR)[2],dim(AR)[3], 1),
                              level = compression
                              )
              h5write(AR, file=paste0(i, ".h5"),
                      name="stacked_channels")
              h5ls(paste0(i, ".h5"))
            }


        ## Plots
            
            if(isTRUE(plots)){
              
              message("  -- creating plots")
              
              if(!is.null(random.crop.x)){
                if(!is.null(random.crop.y)){
                  
                  dir.create('HDF5 cropped plots')
                  dir.create(paste0('HDF5 cropped plots/', i))
                  
                  for(a in channels){
                    make.spatial.plot(roi.dat,
                                      image.roi = i,
                                      image.channel = a,
                                      image.min.threshold = 0,
                                      image.max.threshold = 1,
                                      image.y.flip = FALSE, 
                                      path = paste0('HDF5 cropped plots/', i)
                    )
                  }
                }
                
              } else {
                
                dir.create('HDF5 plots')
                dir.create(paste0('HDF5 plots/', i))
                
                for(a in channels){
                  make.spatial.plot(roi.dat,
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

            
            ## Message
            if(!is.null(random.crop.x)){
              if(!is.null(random.crop.y)){
                
                message('  -- cropped HDF5 file for ', i, ' complete')
                message(' ')
                
              }
            } else {
              message('  -- HDF5 file for ', i, ' complete')
              message(' ')
            }
            
            rm(roi.dat)
            rm(AR)
            gc()
      }
      
      gc()
      
}

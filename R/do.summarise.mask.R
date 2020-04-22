#' do.summarise.mask
#'
#' @export

do.summarise.mask <- function(dat, # list of rasters of one ROI
                              polygon, # polygon (mask) raster
                              fun = "mean"
){

  message("This is a developmental Spectre-spatial function that is still in testing phase with limited documentation. We recommend only using this function if you know what you are doing.")

  ## Create a list of 'raster' names for measurement
      #dat <- spatial.list[["TAXXX 20190705 Liver 243788_s1_p1_r1_a1_ac"]]
      raster.names <- names(dat)
      #raster.names <- raster.names[c(1:3)]

  ## Setup initial mask summary data (X & Y centroids coordinates)
      #polygon <- mask_01_outline
      ply.df <- as.data.frame(polygon)
      ply.df

      ply.centroids <- gCentroid(polygon,byid=TRUE)
      ply.centroids.df <- as.data.frame(ply.centroids)
      ply.centroids.df # mask number, with X and Y coordinates

  ## Loop to measure information within centroid (can be sped up by parallelising? Or maybe better to just do these measurements in CellProfiler?)

      message("Measuring expression information within cell masks. This might take some time, please be patient")

      for(i in raster.names){
        message(paste0("Calculating expression data within cell masks for ", i))
        # i <- raster.names[[1]]
        temp.dat <- dat[[i]]

        extracted.dat = raster::extract(x = temp.dat, y = polygon, df = TRUE) # this is the time consuming step
        extracted.dat.res <- aggregate(. ~ID, data = extracted.dat, FUN = fun)
        colnames(extracted.dat.res)[2] <- i # should we be removing .tiff here? If we do should be the same in the other read.spatial function, to ensure matching consistency

        ply.centroids.df <- cbind(ply.centroids.df, extracted.dat.res[2]) ## doing this would remove the necessity to calculate centroids within the 'make.spatial.plot' function
      }

  return(ply.centroids.df)
}




# https://geocompr.robinlovelace.net/geometric-operations.html

    # rst <- spatial.list[["TAXXX 20190705 Liver 243788_s1_p1_r1_a1_ac"]][["168Er-CD8a-titration_Er168.tiff"]]
    #
    #
    # test = raster::extract(x = rst, y = ply, df = TRUE)
    # test
    #
    # test.res <- aggregate(. ~ID, data = test, FUN = "mean")
    # colnames(test.res)[2] <- "DNA_01"
    # test.res
    #
    # ###
    #
    # ttl <- cbind(ply.centroids.df, test.res[2]) ## doing this would remove the necessity to calculate centroids within the 'make.spatial.plot' function
    # ttl




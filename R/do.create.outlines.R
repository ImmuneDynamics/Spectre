#' do.create.outlines
#'
#' @param spatial.dat NO DEFAULT. Spatial data list
#' @param mask.name NO DEFAULT. Name of the mask to create outlines for
#' @param method DEFAULT = 'stars'. Can be 'stars' or 'raster'
#'
#' @import data.table
#'
#' @export

do.create.outlines <- function(spatial.dat,
                               mask.name,
                               method = 'stars' # 'stars' 'raster'
){

  ### Setup

      #message("This is a developmental Spectre-spatial function that is still in testing phase with limited documentation. We recommend only using this function if you know what you are doing.")

      require(raster)
      require(tiff)
      require(rgeos)
      require(tidyr)
      require(ggplot2)
      require(dplyr)

      # polygons.name <- paste0(mask.name, "_polygons")
      # outlines.name <- paste0(mask.name, "_outlines")
      # centroids.name <- paste0(mask.name, "_centroids")

  ### Initial warning
  
      os.deets <- session_info()
      if(grepl('Windows', os.deets$platform$os)){
        message("Warning: the generation of mask outlines may differ between Mac and Windows. See 'https://github.com/ImmuneDynamics/Spectre/issues/56' for more information.")
      }
      
  ### Slow or fast version

      if(method == 'stars'){
        message(paste0("Creating polygons, outlines, and centroids using 'stars' method."))
        require(stars)
        require(sf)
        require(sp)
        require(s2)
      }

      if(method == 'raster'){
        message(paste0("Creating polygons, outlines, and centroids using standard method -- this step may take some time, please be patient"))
      }

      if(method == 'gdal'){
        message("GDAL version not currently supportedb -- reverting to stars method")
        # if(length(Sys.which("gdal_polygonize.py")) > 1){
        #   message(paste0("Creating polygons, outlines, and centroids using GDAL -- this step may take some time, please be patient"))
        # }
      }

  ### Run

      for(i in names(spatial.dat)){
        # i <- names(spatial.dat)[[1]]
        start.time <- Sys.time()

        mask <- spatial.dat[[i]]$MASKS[[mask.name]]$maskraster

        message(paste0("  -- Processing masks for ROI ", i))

        ## rasterToPolygons method
            if(method == 'raster'){
              polygon <- rasterToPolygons(mask, dissolve=TRUE) # This is the long step
              spatial.dat[[i]]$MASKS[[mask.name]][["polygons"]] <- polygon
              message("    ... polygons complete")
            }

        ## stars method
            if(method == 'stars'){

              require(stars)
              require(sf)
              require(sp)
              require(s2)

              names(mask) <- "TEMP_MASK"

              stars.mask <- stars::st_as_stars(mask)

              #sf::st_crs(stars.mask) <- 4326

              res <- sf::st_as_sf(stars.mask, # requires the sf, sp, raster and stars packages
                                  as_points = FALSE,
                                  merge = TRUE) #,
                                  
                                  #na.rm = TRUE)
                                  #group = TRUE) # TRUE crashes, FALSE does not

              res$TEMP_MASK

              res <- st_make_valid(res)
              
              res <- res %>%
                group_by(TEMP_MASK) %>%
                summarise(geometry = sf::st_union(geometry)) %>%
                ungroup()

              polygon <- sf::as_Spatial(res)

              names(polygon) <- mask.name
              crs(polygon) <- NA

              spatial.dat[[i]]$MASKS[[mask.name]][["polygons"]] <- polygon
              message("     ... polygons complete")
            }

    ### Create outlines
        suppressMessages(
          outline <- fortify(polygon)
        )
        spatial.dat[[i]]$MASKS[[mask.name]][["outlines"]] <- outline
        message("     ... outlines complete")

    ### Create centroids
        centroids <- gCentroid(polygon,byid=TRUE)
        spatial.dat[[i]]$MASKS[[mask.name]][["centroids"]] <- centroids
        message("     ... centroids complete")
      }

  message("Returning spatial data")
  return(spatial.dat)
}

#' read.spatial.files - read in TIFF files into a spatial data object (list of rasters)
#'
#' @export

read.spatial.files <- function(roi.loc, # Corresponding directory of ROI folders
                               rois,

                               mask.loc = NULL,
                               masks = NULL,
                               mask.ext = NULL,

                               cell.dat = NULL, # data.table
                               image.dat = NULL, # data.table
                               panel.dat = NULL, # data.table

                               ext = ".tiff",
                               convert.to.raster = TRUE,
                               correct.extent = TRUE,
                               value.modifier = 65535,
                               calculate.polygons = TRUE,
                               calculate.outlines = TRUE,
                               calculate.centroids = TRUE)
{

    ### Testing data

          # library(Spectre)
          # package.load()
          #
          # library(tiff)
          # library(raster)
          # library(rgeos)
          #
          # ext = ".tiff"
          # convert.to.raster = TRUE
          # correct.extent = TRUE
          # value.modifier = 65535
          # calculate.polygons = TRUE
          # calculate.outlines = TRUE
          # calculate.centroids = TRUE
          #
          # setwd("/Users/thomasa/Desktop/Desktop_items/Modified_Bodenmiller_pipeline/output/histocat/")
          #
          # roi.loc = getwd()
          # rois <- list.dirs(roi.loc, full.names = FALSE, recursive = FALSE)
          # rois <- c(rois[c(6,14)])
          #
          # setwd("/Users/thomasa/Desktop/Desktop_items/Modified_Bodenmiller_pipeline/measurement/")
          #
          # mask.loc <- getwd()
          # masks <- list.files(mask.loc, pattern = ext)
          # mask.ext <- "_ilastik_s2_Probabilities_mask.tiff"
          #
          # setwd("/Users/thomasa/Desktop/Desktop_items/Modified_Bodenmiller_pipeline/1_panel/")
          # panel.dat <- fread("panel.csv")
          #
          # setwd("/Users/thomasa/Desktop/Desktop_items/Modified_Bodenmiller_pipeline/measurement/")
          #
          # cell.dat <- fread("cell.csv")
          # image.dat <- fread("Image.csv")

    message("This is a developmental Spectre-spatial function that is still in testing phase with limited documentation. We recommend only using this function if you know what you are doing.")

    ### Sort rois
    message("Listing ROIs")

          roi.names <- sort(rois, decreasing = FALSE)
          roi.names

    ### Sort mask names, and make some checks to ensure the names are compatible with the rois
    message("Listing and checking masks")

          masks <- sort(masks, decreasing = FALSE)

          mask.roi.names <- gsub(mask.ext, "", masks)
          mask.roi.names <- sort(mask.roi.names, decreasing = FALSE)
          mask.roi.names

          for(i in c(1:length(roi.names))){
            if(roi.names[[i]] != mask.roi.names[[i]]){
              message("Error: in the read.spatial.files function, we found a mismatch between the names of the ROI folders and the names of the mask files. Please check the names of the ROIs and masks, as well as mask.ext")

              message("ROI list:")
              print(as.matrix(roi.names))

              message("Mask list:")
              print(as.matrix(mask.roi.names))

              stop("Alternatively, you can try reading in the tiff files without reading in masks.")
            }
          }

    ### Read TIFFs into a raster.stack, and add to the ROI list
    message("Reading TIFF files into raster.stack")

        setwd(roi.loc)
        ROI.list <- list()
        tiffs.per.roi <- list()

        for(i in roi.names){

          ## Setup
              #i <- roi.names[1]
              setwd(roi.loc)
              setwd(i)

              # tiffs.per.roi[[i]] <- list.files(getwd(), ext)
              # if(!exists("tiffs.per.roi[[i]]")){
              #   stop(paste0("We could not find any ", ext, " files in the first ROI directory (", i, ")"))
              # }

              TIFFs <- list.files(getwd(), ext)
              TIFF.list <- list()

          ## Loop to read in each tiff
              for(a in TIFFs){
                #a <- TIFFs[1]
                TIFF.list[[a]] <- readTIFF(a)

                if(convert.to.raster == TRUE){
                  tiff.map <- TIFF.list[[a]]
                  tiff.raster <- raster(tiff.map)

                  if(correct.extent == TRUE){
                    extent(tiff.raster) <- c(0, dim(tiff.raster)[2], 0,dim(tiff.raster)[1]) # Y axis - X axis
                  }

                  values(tiff.raster) <- values(tiff.raster)*value.modifier
                  tiff.raster <- flip(tiff.raster, 'y')
                  TIFF.list[[a]] <- tiff.raster
                }
              }

          ## Remove mask from list of rasters, if it is present
              # val <- grep(mask.ext, names(TIFF.list))
              # mask <- TIFF.list[val]
              # TIFF.list <- TIFF.list[-val]

              for(n in c(1:length(names(TIFF.list)))){
                names(TIFF.list)[n] <- gsub(ext, "", names(TIFF.list)[n])
              }

          ## Finalise for ROI

              raster.stack <- stack(TIFF.list)
              #ROI.list[[i]][["tiffs"]] <- raster.stack
              ROI.list[[i]]$rasters <- raster.stack
              ROI.list[[i]]$raster.names <- names(raster.stack)
        }

        # ROI.list$`20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac`$rasters
        # ROI.list$`20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac`$raster.names
        #

    ### Read mask TIFFs into a list, then add the masks to each ROI
    message("Reading mask tiffs as a raster")

        setwd(mask.loc)
        mask.list <- list()
        mask.list

        for(i in masks){
          # i <- masks[1]
          mask.list[[i]] <- readTIFF(i)
          mask.list[[i]] <- raster(mask.list[[i]])

          if(correct.extent == TRUE){
            extent(mask.list[[i]]) <- c(0, dim(mask.list[[i]])[2], 0,dim(mask.list[[i]])[1]) # Y axis - X axis
          }
          mask.list[[i]] <- flip(mask.list[[i]], 'y')
        }

        names(mask.list) <- gsub(mask.ext, "", names(mask.list))
        names(mask.list)

        for(i in c(1:length(names(ROI.list)))){
          # i <- c(1:length(names(ROI.list)))[1]
          a <- names(ROI.list)[i]
          b <- mask.roi.names[[i]]

          if(a != b){
            stop("The ROI name for the rasters and the mask do not match")
          }

          ROI.list[[i]]$masks <- mask.list[[i]]
        }

    ### Review

        # ROI.list$`20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac`$rasters
        # ROI.list$`20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac`$raster.names
        # ROI.list$`20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac`$masks

    ### Calculate polygons, outlines, and centroid coordinates (might take a few minutes)
    message("Calculating mask polygons, outlines, and centroids")

      for(i in names(ROI.list)){
        # i <- names(ROI.list)[[1]]

        if(calculate.polygons == TRUE){
          message(paste0("Calculating mask polygons for ", i))
          ROI.list[[i]]$polygons <- do.convert.mask.to.outline(ROI.list[[i]]$masks)
        }

        if(calculate.outlines == TRUE){
          outline <- fortify(ROI.list[[i]]$polygons)
          ROI.list[[i]]$outlines <- outline
        }

        if(calculate.centroids == TRUE){
          temp <- gCentroid(ROI.list[[i]]$polygons,byid=TRUE)
          #temp <- as.data.frame(temp)
          ROI.list[[i]]$centroids <- temp
        }

      }

    ### Review

        # ROI.list$`20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac`$rasters
        # ROI.list$`20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac`$raster.names
        # ROI.list$`20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac`$masks
        # ROI.list$`20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac`$polygons
        # ROI.list$`20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac`$outlines
        # ROI.list$`20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac`$centroids

    ### Finalise
    message("Finalising ROI and mask addition to spatial data")

        spatial.dat <- list()
        spatial.dat$rois <- ROI.list

        spatial.dat$meta.data$roi.names <- names(ROI.list)

        all.raster.names <- list()

        for(i in names(ROI.list)){
          all.raster.names[[i]] <- ROI.list[[i]]$raster.names
        }

        spatial.dat$meta.data$raster.names <- all.raster.names

        spatial.dat$rois
        spatial.dat$meta.data

    ### Add in cell and image data
    message("Adding 'per cell' data extracted in CellProfiler")

        ## CP outputs

            measure.dat <- cell.dat
            as.matrix(names(measure.dat))
            dim(measure.dat)

            spatial.dat$cp.cell <- measure.dat

            image.dat <- image.dat
            as.matrix(names(image.dat))
            image.dat

            spatial.dat$cp.image <- image.dat

        ## Panel

            spatial.dat$panel <- panel.dat
            panel <- panel.dat

            new.names <- panel[panel[["full"]] == 1,]$Target
            new.names

        ##
            measure.dat$ImageNumber
            image.dat$FileName_CellImage

            measure.dat <- do.embed.columns(measure.dat,
                                            type = "data.table",
                                            base.name = "ImageNumber",
                                            col.name = "ImageName",
                                            match.to = unique(image.dat$ImageNumber),
                                            new.cols = unique(image.dat$FileName_CellImage))

            measure.dat$ImageName
            names.length <- length(names(measure.dat))
            names.length
            names <- names(measure.dat)[c(1:(names.length-1))]
            names

            measure.dat <- data.table("ImageName" = measure.dat$ImageName, measure.dat[,names,with = FALSE])
            as.matrix(names(measure.dat))
            dim(measure.dat)


        ### Filter data to most relevant cell data

            as.matrix(names(measure.dat))

            nms <- measure.dat[,c(1:3),]
            nms <- data.table(nms, "X" = measure.dat$Location_Center_X, "Y" = measure.dat$Location_Center_Y)
            nms

        ### Extract mean intensity

            means <- measure.dat[,grepl( "Intensity_MeanIntensity_FullStack_" , names(measure.dat)),with = FALSE]
            means <- means*65535
            names(means) <- gsub("Intensity_MeanIntensity_FullStack_c", "", names(means))
            names(means)

            for(i in c(1:length(names(means)))){
              a <- names(means)[i]
              if(nchar(a) == 1){
                a <- paste0("0", a)
                names(means)[i] <- a
              }
            }

            means <- means[,sort(names(means)),with = FALSE]
            names(means)
            names(means) <- new.names


        ### Extract mean intensity (filtered)

            means.filtered <- measure.dat[,grepl( "Intensity_MeanIntensity_FullStackFiltered_" , names(measure.dat)),with = FALSE]
            means.filtered <- means.filtered*65535
            names(means.filtered) <- gsub("Intensity_MeanIntensity_FullStackFiltered_c", "", names(means.filtered))
            names(means.filtered)

            for(i in c(1:length(names(means.filtered)))){
              a <- names(means.filtered)[i]
              if(nchar(a) == 1){
                a <- paste0("0", a)
                names(means.filtered)[i] <- a
              }
            }

            means.filtered <- means.filtered[,sort(names(means.filtered)),with = FALSE]
            names(means.filtered)
            names(means.filtered) <- new.names

        ### Fix up ROI name (take off the mask.ext)

            for(i in unique(nms$ImageName)){
              new.name <- gsub(mask.ext, "", i)

              temp <- nms[nms[["ImageName"]] == i,]
              nrow(temp)

              new.names.vector <- rep(new.name, nrow(temp))

              nms[nms[["ImageName"]] == i,]$ImageName <- new.names.vector
            }

        ### Add to spatial.dat

            for(i in unique(nms$ImageName)){
              # i <- "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac"
              means.temp <- cbind(nms, means)
              spatial.dat$rois[[i]]$per.cell.means <- means.temp[means.temp[["ImageName"]] == i,]

              means.filtered.temp <- cbind(nms, means.filtered)
              spatial.dat$rois[[i]]$per.cell.means.filtered <- means.filtered.temp[means.filtered.temp[["ImageName"]] == i,]
            }

    ### Process complete
        message("Process COMPLETE")

  return(spatial.dat)
}


#' make.spatial.plot
#' 
#' @param dat NO DEFAULT.A spatial data list
#' @param image.roi NO DEFAULT.The name of the ROI to plot
#' @param image.channel NO DEFAULT.Name of the channel to plot
#' @param mask.outlines DEFAULT = NULL. Name of the mask for outlines
#' @param cell.dat DEFAULT = NULL. Character name for the cellular dataset in the spatial data list, or a data.table object.
#' @param cell.col DEFAULT = NULL. Name of cellular dataset to colour points by
#' @param image.y.flip DEFAULT = TRUE. Flips the Y-axis orientation.
#' @param image.mask.size DEFAULT = 0.1
#' @param image.mask.colour DEFAULT = "gold"
#' @param image.min.threshold DEFAULT = 0.00. Lower threshold for plot signal. Values below this target will be clipped.
#' @param image.max.threshold DEFAULT = 0.99. Upper threshold for plot signal. Values above this target will be clipped.
#' @param image.blank DEFAULT = FALSE. Will blank the plot image.
#' @param cell.x DEFAULT = "x". Column name for the 'x' coordinate in data.table
#' @param cell.y DEFAULT = "y". Column name for the 'y' coordinate in data.table
#' @param cell.col.type DEFAULT= "numeric". Can be 'factor'.
#' @param cell.colours = DEFAULT = "spectral".
#' @param cell.col.min.threshold = DEFAULT = 0.01. Values below this target will be clipped.
#' @param cell.col.max.threshold = DEFAULT = 0.995. Values above this target will be clipped.
#' @param title DEFAULT = paste0(image.roi)
#' @param dot.size DEFAULT = 1
#' @param dot.alpha DEFAULT = 1
#' @param align.xy.by DEFAULT = cell.dat
#' @param align.col.by DEFAULT = cell.dat
#' @param save.to.disk DEFAULT = TRUE
#' @param path DEFAULT = getwd(). Path for saving image files
#' @param plot.width DEFAULT = 9
#' @param plot.height DEFAULT = 7
#' @param blank.axis DEFAULT = FALSE
#'
#' @import data.table
#'
#' @export

make.spatial.plot <- function(dat, # spatial data object
                              image.roi, # name of ROI
                              image.channel, # name of channel

                              ## Options for adding cell outlines
                              mask.outlines = NULL, # character -- the outlines in dat object

                              ## Options for adding cellular data
                              cell.dat = NULL, # can be character (if it's data within dat) or a data.table
                              cell.col = NULL, # column for colouration

                              ## Other settings (with defaults)
                              image.y.flip = TRUE,
                              image.mask.size = 0.1,
                              image.mask.colour = "gold",
                              image.min.threshold = 0.00,
                              image.max.threshold = 0.99,
                              image.blank = FALSE,

                              cell.x = "x",
                              cell.y = "y",
                              cell.col.type = "numeric",
                              cell.colours = "spectral",
                              cell.col.min.threshold = 0.01,
                              cell.col.max.threshold = 0.995,

                              title = paste0(image.roi),
                              dot.size = 1,
                              dot.alpha = 1,
                              align.xy.by = cell.dat, # choose a data frame to set absolute limits for X/Y/colour
                              align.col.by = cell.dat,
                              save.to.disk = TRUE,
                              path = getwd(),
                              plot.width = 9,
                              plot.height = 7,
                              blank.axis = FALSE)
{

  ### TESTING
      # library(raster)
      # library(data.table)
      # library(tiff)
      # library(ggplot2)
      #
      # dat = dat
      #
      # dat$meta.data
      #
      # roi = "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac"
      # roi.marker = "CD20_Dy161"
    
      # cell.dat <- dat$cell.dat.means.filtered
      # cell.dat <- cell.dat[cell.dat[["ImageName"]] == "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac_ilastik_s2_Probabilities_mask.tiff",]
      # cell.dat = cell.dat
      # cell.x = "X"
      # cell.y = "Y"
      # cell.colour = 'CD20'
      # 
      # add.outlines = TRUE
      # flip.y.axis = TRUE

  ### Check that necessary packages are installed
      if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
      if(!is.element('ggplot2', installed.packages()[,1])) stop('ggplot2 is required but not installed')
      if(!is.element('scales', installed.packages()[,1])) stop('scales is required but not installed')
      if(!is.element('colorRamps', installed.packages()[,1])) stop('colorRamps is required but not installed')
      if(!is.element('ggthemes', installed.packages()[,1])) stop('ggthemes is required but not installed')
      if(!is.element('RColorBrewer', installed.packages()[,1])) stop('RColorBrewer is required but not installed')
      if(!is.element('raster', installed.packages()[,1])) stop('raster is required but not installed')
      
  ### Require packages
      require(ggplot2)
      require(scales)
      require(colorRamps)
      require(ggthemes)
      require(RColorBrewer)
      require(raster)
      
  ### Compatability conversions
    
      roi <- image.roi
      roi.marker <- image.channel
    
      #cell.dat
      # if('data.table' %in% class(cell.dat)){
      #   cell.dat <- cell.dat[cell.dat[['ROI']] == image.roi,]
      # }
      
      cell.colour <- cell.col
    
      add.outlines <- image.outlines <- mask.outlines
      flip.y.axis <- image.y.flip
    
      cell.colour.type <- cell.col.type
    
      raster.mask.size <- image.mask.size
      raster.mask.colour <- image.mask.colour
      raster.min.threshold <- image.min.threshold
      raster.max.threshold <- image.max.threshold
    
      col.min.threshold <- cell.col.min.threshold
      col.max.threshold <- cell.col.max.threshold
    
      colours <- cell.colours

  ### Colour setup
    
      # Jet
      if(colours == "jet"){
        colour.scheme <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
      }
    
      # Spectral
      if(colours == "spectral"){
        spectral.list <- colorRampPalette(brewer.pal(11,"Spectral"))(50)
        spectral.list <- rev(spectral.list)
        colour.scheme <- colorRampPalette(c(spectral.list))
      }
    
      # Viridis
      if(colours == "viridis"){
        colour.scheme <- colorRampPalette(c(viridis_pal(option = "viridis")(50)))
      }
    
      # Inferno
      if(colours == "inferno"){
        colour.scheme <- colorRampPalette(c(viridis_pal(option = "inferno")(50)))
      }
    
      #Magma
      if(colours == "magma"){
        colour.scheme <- colorRampPalette(c(viridis_pal(option = "magma")(50)))
      }

  ### cell.dat setup
    
      if(!is.null(cell.dat)){
    
        if(is.character(cell.dat) == TRUE){
          temp <- dat[[roi]]@DATA[[cell.dat]]
          cell.dat <- temp
        }
    
    
        if(cell.colour.type == "numeric"){
    
          # Dot point colouration
          if(is.null(align.col.by) == TRUE){
            ColrMin <- quantile(cell.dat[[cell.colour]], probs = c(col.min.threshold))
            ColrMax <- quantile(cell.dat[[cell.colour]], probs = c(col.max.threshold))
          }
    
          if(is.null(align.col.by) == FALSE){
            ColrMin <- quantile(align.col.by[[cell.colour]], probs = c(col.min.threshold))
            ColrMax <- quantile(align.col.by[[cell.colour]], probs = c(col.max.threshold))
          }
    
        }
    
      }


  ### Preparat the raster data
    
      ## Image prep
    
      raster.image <- dat[[roi]]@RASTERS[[roi.marker]]
    
      tiff.p <- rasterToPoints(raster.image)
      tiff.df <- data.frame(tiff.p)
      raster.label <- names(tiff.df)[3]
      colnames(tiff.df) <- c("x_axis", "y_axis", raster.label)
    
      ## Create cell outlines
      if(!is.null(mask.outlines)){
        outline <- dat[[roi]]@MASKS[[mask.outlines]]$outlines
        centroids <- dat[[roi]]@MASKS[[mask.outlines]]$centroids
    
        centroid.xmin <- min(centroids[,1])
        centroid.xmax <- max(centroids[,1])
    
        centroid.ymin <- min(centroids[,2])
        centroid.ymax <- max(centroids[,2])
      }
    
      ## Flip y-axis values
    
      # if(flip.y.axis == TRUE){
      #   dat <- invert.y.axis(dat, y.axis)
      # }
    
      ## Normalise XY for cell centroids
      plot.normalize <- function(dat, min, max){
        return(((dat- min(dat)) / (max(dat)-min(dat))) * (max - min) + min)
      }
    
      # if(!is.null(cell.dat)){
      #   # X AXIS
      #   cell.dat[[cell.x]] <- plot.normalize(cell.dat[[cell.x]], min = centroid.xmin, max = centroid.xmax)
      #
      #   # Y AXIS
      #   cell.dat[[cell.y]] <- plot.normalize(cell.dat[[cell.y]], min = centroid.ymin, max = centroid.ymax)
      #
      # }
    
      ## Raster colour limits
    
      RastMin <- quantile(tiff.df[[3]], probs = c(raster.min.threshold))
      RastMax <- quantile(tiff.df[[3]], probs = c(raster.max.threshold))

  ###############################################
  ### Add a check to see if centroids line up ###
  ###############################################

  ### Generate and show coloured plot
    
      if(image.blank == FALSE){
        p <- ggplot(data=tiff.df, aes(x=tiff.df[[1]], y=tiff.df[[2]])) +
    
          ## Plot the raster (IMC image)
          geom_raster(aes(fill=tiff.df[[3]])) +
          scale_fill_gradient(raster.label,
                              low = "black",
                              high = "white",
                              limits=c(RastMin,RastMax),
                              oob=squish)
      }
    
      if(image.blank == TRUE){
        p <- ggplot(data=tiff.df, aes(x=tiff.df[[1]], y=tiff.df[[2]])) +
    
          ## Plot the raster (IMC image)
          geom_raster(aes(fill=tiff.df[[3]])) +
          scale_fill_gradient(raster.label,
                              low = "black",
                              high = "black",
                              limits=c(RastMin,RastMax),
                              oob=squish)
      }
    


  ### Plot the cell mask boundaries
    
      if(!is.null(image.outlines)){
        p <- p + geom_path(aes(x = long, y = lat, group = group),
                           data = outline,
                           size = raster.mask.size,
                           col = raster.mask.colour)
      }
    
      ## Plot the cellular data
    
      if(!is.null(cell.dat)){
        if(cell.colour.type == "numeric"){
          p <- p + geom_point(data=cell.dat,
                              aes(x=cell.dat[[cell.x]], y=cell.dat[[cell.y]], color = cell.dat[[cell.colour]]),  #as.numeric(as.character(col))
                              size = dot.size, #dot.size
                              alpha = dot.alpha # shape = 1
          ) +
    
            scale_color_gradientn(colours = colour.scheme(50),
                                  limits = (c(ColrMin,ColrMax)),
                                  oob=squish,
                                  name = cell.colour)
        }
    
        if(cell.colour.type != "numeric"){
          p <- p + geom_point(data=cell.dat,
                              aes(x=cell.dat[[cell.x]], y=cell.dat[[cell.y]], color = as.factor(cell.dat[[cell.colour]])),  #as.numeric(as.character(col))
                              size = dot.size, #dot.size
                              alpha = dot.alpha # shape = 1
          ) +
    
            scale_colour_discrete(name = cell.colour)
        }
    
      }
    
      ## Setup some themes
      p <- p + theme_bw() +
        coord_equal() +
        xlab(cell.x)+
        ylab(cell.y)+
        ggtitle(title)
    
      ## More themes
      p <- p + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # change 'colour' to black for informative axis
                     axis.title.x=element_text(color="Black", face="bold", size=18),
                     axis.title.y=element_text(color="Black", face="bold", size=18),
                     legend.text=element_text(size=12), # large = 30 # small = 8
                     legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
                     legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
                     #legend.title=element_blank(),
                     plot.title = element_text(color="Black", face="bold", size=16, hjust=0) # size 70 for large, # 18 for small
      )
    
      if(flip.y.axis == TRUE){
    
        p <- p + scale_y_reverse()
    
      }

  ### Save ggplot to disk if desired
      if(save.to.disk == TRUE){
        ggsave(filename = paste0(title, "_ROI_", roi.marker, "_marker_", cell.colour,".png"),
               plot = p,
               path = path,
               width = plot.width,
               height = plot.height,
               limitsize = FALSE)
      }
    
      print(p)

}






#' make.colour.plot - Create a dot plot (X vs Y) coloured by a selected continuous column (e.g. marker expression)
#'
#' This function allows you to create a coloured XY plot where each cell is coloured by a selected column. Typically used to plot cells on tSNE1/2 or UMAP1/2 coloured by select cellular markers.
#'
#' @usage make.colour.plot(d, x, y, col.axis)
#'
#' @param dat NO DEFAULT. data.table Input sample.
#' @param x.axis NO DEFAULT. Character. Column for X axis.
#' @param y.axis NO DEFAULT. Character. Column for Y axis.
#' @param col.axis NO DEFAULT. Character. Column for colour.
#' @param title DEFAULT = col.axis. Character. Title for the plot.
#' @param col.min.threshold DEFAULT = 0.01. Numeric. Define minimum threshold for colour scale. Values below this limit will be coloured as the chosen minimum threshold.
#' @param col.max.threshold DEFAULT = 0.995 Numeric. Define maximum threshold for colour scale. Values above this limit will be coloured as the chosen maximum threshold.
#' @param colours DEFAULT = "spectral". Character. Colour scheme to use. Can be "spectral", "jet", "viridis", "magma", or "inferno".
#' @param dot.size DEFAULT = 1. Numeric. Size of the dots.
#' @param align.xy.by DEFAULT = dat. data.table Sample to use to determine minimum and maximum X and Y axis values.
#' @param align.col.by DEFAULT = dat. data.table. Sample to use to determine minimum and maximum colour values.
#' @param save.to.disk DEFAULT = TRUE. Will save the ggplot to disk. If FALSE, will only show the ggplot.
#' @param path DEFAULT = getwd(). The location to save your ggplot. By default, will save to current working directory. Can be overidden.
#' @param plot.width DEFAULT = 9. Width of the ggplot when saved to disk.
#' @param plot.height DEFAULT = 7. Height of the ggplot when saved to disk.
#' @param blank.axis DEFAULT = FALSE Logical, do you want a minimalist graph?
#'
#' @return prints and saves a ggplot.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @usage See \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#'
#' @examples
#' make.colour.plot(d = demo.umap, x = "UMAP_X", y = "UMAP_Y", col.axis = "BV605.Ly6C")
#'
#' @export

make.colour.plot <- function(dat,
                             x.axis,
                             y.axis,
                             col.axis,
                             title = col.axis,
                             col.min.threshold = 0.01,
                             col.max.threshold = 0.995,
                             colours = "spectral",
                             dot.size = 1,
                             align.xy.by = dat, # choose a data frame to set absolute limits for X/Y/colour
                             align.col.by = dat,
                             save.to.disk = TRUE,
                             path = getwd(),
                             plot.width = 9,
                             plot.height = 7,
                             blank.axis = FALSE)

{

  ## TESTING
      # dat <- demo.umap
      # x.axis = "UMAP_42_X"
      # y.axis = "UMAP_42_Y"
      # col.axis = "BV605.Ly6C"
      # bubble.lab = "FlowSOM_metacluster"
      # col.min.threshold = 0.01
      # col.max.threshold = 1.0
      # title = paste0("All samples", " - ", "BV605.Ly6C")
      # colours = "spectral"
      # dot.size = 1
      # align.xy.by = dat
      # align.col.by = dat

  ## Check that necessary packages are installed
    if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required by not installed')
    if(!is.element('ggplot2', installed.packages()[,1])) stop('ggplot2 is required by not installed')
    if(!is.element('scales', installed.packages()[,1])) stop('scales is required by not installed')
    if(!is.element('colorRamps', installed.packages()[,1])) stop('colorRamps is required by not installed')
    if(!is.element('ggthemes', installed.packages()[,1])) stop('ggthemes is required by not installed')
    if(!is.element('RColorBrewer', installed.packages()[,1])) stop('RColorBrewer is required by not installed')

  ## Require packages
    require(Spectre)
    require(ggplot2)
    require(scales)
    require(colorRamps)
    require(ggthemes)
    require(RColorBrewer)

  ## Colour setup

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

  ## Define limits

      # X AXIS
      if(is.null(align.xy.by) == TRUE){
        Xmax <- max(dat[[x.axis]])
        Xmin <- min(dat[[x.axis]])
      }

      if(is.null(align.xy.by) == FALSE){
        Xmax <- max(align.xy.by[[x.axis]])
        Xmin <- min(align.xy.by[[x.axis]])
        }


      # Y AXIS
      if(is.null(align.xy.by) == TRUE){
        Ymax <- max(dat[[y.axis]])
        Ymin <- min(dat[[y.axis]])
        }

      if(is.null(align.xy.by) == FALSE){
        Ymax <- max(align.xy.by[[y.axis]])
        Ymin <- min(align.xy.by[[y.axis]])
        }

      # COLOUR
      if(is.null(align.col.by) == TRUE){
        ColrMin <- quantile(dat[[col.axis]], probs = c(col.min.threshold))
        ColrMax <- quantile(dat[[col.axis]], probs = c(col.max.threshold))
      }

      if(is.null(align.col.by) == FALSE){
        ColrMin <- quantile(align.col.by[[col.axis]], probs = c(col.min.threshold))
        ColrMax <- quantile(align.col.by[[col.axis]], probs = c(col.max.threshold))
      }


  ## Generate and show coloured plot

    p <- ggplot(data = dat, aes(x = dat[[x.axis]], y = dat[[y.axis]], colour = dat[[col.axis]])) +
                geom_point(size = dot.size) +

                scale_colour_gradientn(colours = colour.scheme(50),
                                       limits = c(ColrMin, ColrMax),
                                       oob=squish) +
                ggtitle(title) +
                xlim(Xmin, Xmax) +
                ylim(Ymin, Ymax) +

                xlab(x.axis)+
                ylab(y.axis)

  ## Add some themes

    p <- p + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # change 'colour' to black for informative axis
                   axis.title.x=element_text(color="Black", face="bold", size=18),
                   axis.title.y=element_text(color="Black", face="bold", size=18),
                   legend.text=element_text(size=12), # large = 30 # small = 8
                   legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
                   legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
                   legend.title=element_blank(),
                   plot.title = element_text(color="Black", face="bold", size=22, hjust=0) # size 70 for large, # 18 for small
                   )

  ## Blank the axes if desired

    if(blank.axis == TRUE){
      p <- p + theme(axis.line=element_blank(),
                     axis.text.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     panel.grid.major = element_blank(),
                     panel.background=element_blank(),
                     panel.border=element_blank(),
                     panel.grid.minor=element_blank(),
                     plot.background=element_blank(),
                     # legend.position = "right",
                     # legend.text=element_text(size=15), # large = 30 # small = 8
                     # legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
                     # legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
                     #legend.title=element_blank(),
                     #plot.title = element_text(color="Black", face="bold", size=15, hjust=0
                     )
      }

  ## Save ggplot to disk if desired
      if(save.to.disk == TRUE){
        ggsave(filename = paste0(title, ".png"),
               plot = p,
               path = path,
               width = plot.width,
               height = plot.height,
               limitsize = FALSE)
      }

      print(p)

}

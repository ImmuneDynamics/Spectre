#' make.density.plot - ...
#'
#' @usage make.density.plot(dat, ...)
#'
#' @param dat NO DEFAULT. data.frame.
#' @param x.axis NO DEFAULT. Character. Column for X axis.
#' @param y.axis NO DEFAULT. Character. Column for Y axis.
#' @param title DEFAULT = "Density". Character. Title for the plot.
#' @param colours DEFAULT = "spectral". Character. Colour scheme to use. Can be "spectral", "jet", "viridis", "magma", or "inferno".
#' @param dot.size DEFAULT = 1. Numeric. Size of the dots.
#' @param align.xy.by DEFAULT = dat. data.table Sample to use to determine minimum and maximum X and Y axis values.
#' @param align.col.by DEFAULT = dat. data.table. Sample to use to determine minimum and maximum colour values.
#' @param save.to.disk DEFAULT = TRUE. Will save the ggplot to disk. If FALSE, will only show the ggplot.
#' @param path DEFAULT = getwd(). The location to save your ggplot. By default, will save to current working directory. Can be overidden.
#' @param plot.width DEFAULT = 9. Width of the ggplot when saved to disk.
#' @param plot.height DEFAULT = 7. Height of the ggplot when saved to disk.
#' @param blank.axis DEFAULT = FALSE. Logical, do you want a minimalist graph?
#'
#' @export

make.density.plot <- function(dat,
                        x.axis,
                        y.axis,
                        title = "Density",
                        colours = "viridis",
                        dot.size = 1,
                        align.xy.by = dat, # choose a data frame to set absolute limits for X/Y/colour
                        align.col.by = dat,
                        save.to.disk = TRUE,
                        path = getwd(),
                        plot.width = 9,
                        plot.height = 7,
                        blank.axis = FALSE
){
  
  ### TEST DATA
  # dat <- Spectre::demo.umap
  # x.axis <- "UMAP_42_X"
  # y.axis <- "UMAP_42_Y"
  
  ## Check that necessary packages are installed
  if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required by not installed')
  if(!is.element('ggplot2', installed.packages()[,1])) stop('ggplot2 is required by not installed')
  if(!is.element('RColorBrewer', installed.packages()[,1])) stop('RColorBrewer is required by not installed')
  if(!is.element('ggpointdensity', installed.packages()[,1])) stop('ggpointdensity is required by not installed')
  
  ## Require packages
  require(Spectre)
  require(ggplot2)
  require(RColorBrewer)
  require(ggpointdensity)
  
  ## Define limits
      # X AXIS
        Xmax <- max(dat[[x.axis]])
        Xmin <- min(dat[[x.axis]])
      
      # Y AXIS
        Ymax <- max(dat[[y.axis]])
        Ymin <- min(dat[[y.axis]])
  
  ## Generate and show coloured plot
      
      p <- ggplot2::ggplot(data = dat, ggplot2::aes(x = dat[[x.axis]], y = dat[[y.axis]])) +
        ggplot2::geom_point(size = dot.size) +
        
        ggplot2::ggtitle(title) +
        ggplot2::xlim(Xmin, Xmax) +
        ggplot2::ylim(Ymin, Ymax) +
        
        ggplot2::xlab(x.axis) +
        ggplot2::ylab(y.axis) +
        
        ggpointdensity::geom_pointdensity()
  
  ## Add some themes
      
      p <- p + ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "black", size = 0.5), # change 'colour' to black for informative axis
                     axis.title.x=ggplot2::element_text(color="Black", face="bold", size=18),
                     axis.title.y=ggplot2::element_text(color="Black", face="bold", size=18),
                     legend.text=ggplot2::element_text(size=12), # large = 30 # small = 8
                     legend.key.height=ggplot2::unit(1,"cm"), # large = 3 # small = 1.2
                     legend.key.width=ggplot2::unit(0.4,"cm"), # large = 1 # small = 0.4
                     legend.title=ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(color="Black", face="bold", size=22, hjust=0) # size 70 for large, # 18 for small
      )
        
        
      if(colours == "viridis" || colours == "magma" || colours == "inferno"){
        p <- p + viridis::scale_colour_viridis(option = colours)
      }
        
      if(colours == "jet") {
        p <- p + ggplot2::scale_colour_gradientn(colours = c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
      }
      
      if(colours == "spectral"){
        p <- p + ggplot2::scale_colour_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(50)))
      }
  
  ## Blank the axes if desired
      
      if(blank.axis == TRUE){
        p <- p + ggplot2::theme(axis.line=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks=ggplot2::element_blank(),
                       axis.title.x=ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.background=ggplot2::element_blank(),
                       panel.border=ggplot2::element_blank(),
                       panel.grid.minor=ggplot2::element_blank(),
                       plot.background=ggplot2::element_blank()
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
        ggplot2::ggsave(filename = paste0(title, ".png"),
               plot = p,
               path = path,
               width = plot.width,
               height = plot.height,
               limitsize = FALSE)
      }
      
      print(p)

}


#ggplot2::ggplot(data = x, mapping = ggplot2::aes(x = x.axis, y = y.axis)) +
#  ggpointdensity::geom_pointdensity() +
#  viridis::scale_colour_viridis()

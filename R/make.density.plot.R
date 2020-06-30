#' Make a density plot (using ggpointdensity::geom_pointdensity())
#'
#' Method to create a density plot.
#' A density plot can summarise the cellular distribution across a given plot.
#' Can be used to plot any two markers.
#' Uses "ggplot2" for plots, "ggpointdensity" for adding density to plots, "RColorBrewer" for colour schemes.
#'
#' @param dat NO DEFAULT. data.frame.
#' @param x.axis NO DEFAULT. Character. Column for X axis.
#' @param y.axis NO DEFAULT. Character. Column for Y axis.

#' @param colours DEFAULT = "viridis". Character. Colour scheme to use. Can be "spectral", "jet", "viridis", "magma", "inferno", or "BuPu".
#' @param dot.size DEFAULT = 1. Numeric. Size of the dots.
#' @param align.xy.by DEFAULT = dat. data.table Sample to use to determine minimum and maximum X and Y axis values.
#' @param title DEFAULT = "Density". Character. Title for the plot.
#' @param dot.size DEFAULT = 1. Numeric. Size of the dots.
#' @param plot.width DEFAULT = 9. Width of the ggplot when saved to disk.
#' @param plot.height DEFAULT = 7. Height of the ggplot when saved to disk.
#' @param nudge_x DEFAULT = 0.5. When add.label = TRUE, distance the label is shifted from the centroid point on the X axis.
#' @param nudge_y DEFAULT = 0.5. When add.label = TRUE, distance the label is shifted from the centroid point on the Y axis.
#' @param square DEFAULT = TRUE. Ensures the plot is saved as a square. Set to FALSE if you want a plot with different X and Y lengths.
#' @param blank.axis DEFAULT = FALSE. Logical, do you want a minimalist graph?
#' @param save.to.disk DEFAULT = TRUE. Will save the ggplot to disk. If FALSE, will only show the ggplot.
#' @param path DEFAULT = getwd(). The location to save your ggplot. By default, will save to current working directory. Can be overidden.
#'
#' @usage make.density.plot(dat, x.axis, y.axis)
#'
#' @examples
# Load packages
#' library(Spectre)
#' package.check()
#' package.load()
#'
#' # Read data
#' cell.dat <- Spectre::demo.umap
#' cell.dat <- as.data.table(cell.dat)
#'
#' # Draw plot
#' Spectre::make.density.plot(dat = cell.dat,
#'                           x.axis = "UMAP_42_X",
#'                           y.axis = "UMAP__42Y")
#'
#' @author Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#'
#' @export

make.density.plot <- function(dat,
                              x.axis,
                              y.axis,

                              colours = "viridis", # can be spectral, jet, etc      # only works for continuous
                              align.xy.by = dat,
                              title = "Density",
                              dot.size = 1,
                              plot.width = 9,
                              plot.height = 7,
                              nudge_x = 0.5,
                              nudge_y = 0.5,
                              square = TRUE,
                              blank.axis = FALSE,
                              save.to.disk = TRUE,
                              path = getwd()
                              ){

  ### Check for packages
  if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
  if(!is.element('ggplot2', installed.packages()[,1])) stop('ggplot2 is required but not installed')
  if(!is.element('RColorBrewer', installed.packages()[,1])) stop('RColorBrewer is required but not installed')
  if(!is.element('ggpointdensity', installed.packages()[,1])) stop('ggpointdensity is required but not installed')

  ### Load packages
  require(Spectre)
  require(ggplot2)
  require(RColorBrewer)
  require(ggpointdensity)

  ### Some tests


  ### Define limits

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

  ### Initialise plot

  p <- ggplot(data = dat,
              aes(x = .data[[x.axis]],
                  y = .data[[y.axis]])) +

    ggpointdensity::geom_pointdensity(size = dot.size) +
    ggtitle(title)

  ### Set up axis
  p <- p + scale_x_continuous(breaks = scales::pretty_breaks(n = 8), name = x.axis, limits = c(Xmin, Xmax))
  p <- p + scale_y_continuous(breaks = scales::pretty_breaks(n = 8), name = y.axis, limits = c(Ymin, Ymax))

  ### Set up themes etc

  p <- p + ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "black", size = 0.5), # change 'colour' to black for informative axis
                          axis.title.x=ggplot2::element_text(color="Black", face="bold", size=18),
                          axis.title.y=ggplot2::element_text(color="Black", face="bold", size=18),
                          legend.text=ggplot2::element_text(size=12), # large = 30 # small = 8
                          legend.key.height=ggplot2::unit(1,"cm"), # large = 3 # small = 1.2
                          legend.key.width=ggplot2::unit(0.4,"cm"), # large = 1 # small = 0.4
                          legend.title=ggplot2::element_blank(),
                          plot.title = ggplot2::element_text(color="Black", face="bold", size=22, hjust=0) # size 70 for large, # 18 for small
  )

  ### Colours

  if(colours == "viridis" || colours == "magma" || colours == "inferno"){
    p <- p + viridis::scale_colour_viridis(option = colours)
  }

  if(colours == "jet") {
    p <- p + ggplot2::scale_colour_gradientn(colours = c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  }

  if(colours == "spectral"){
    p <- p + ggplot2::scale_colour_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(50)))
  }

  #Blue to Purple
  if(colours == "BuPu"){
    colour.list <- (colorRampPalette(RColorBrewer::brewer.pal(9, "BuPu"))(31)) # 256
    #colours <- colorRampPalette(c(colour.list))

    p <- p + ggplot2::scale_colour_gradientn(colours = colour.list)

  }


  ### Blank axis options

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

  if(square == TRUE){
    p <- p+ theme(aspect.ratio=1)
  }

  ### Save plot
  if(save.to.disk == TRUE){

    ggsave(filename = paste0("Density plot - ", title, " - plotted on ", x.axis, " by ", y.axis, ".png"),
           plot = p,
           path = path,
           width = plot.width,
           height = plot.height,
           limitsize = FALSE)
  }

  ### Print plot
  print(p)
}

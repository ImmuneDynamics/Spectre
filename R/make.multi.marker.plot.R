#' Make mulitple plots for each marker
#' 
#' Method to create multiple plots for each marker.
#' This function allows you to create a grid of plots, where the cells are subsetted by a certain factor (e.g. one sample per plot).
#' These can then be coloured by a marker or by another factor (e.g Group).
#' Makes use of Spectre functions and make.colour.plot make.density.plot.
#'
#' @param dat NO DEFAULT. A data frame containing all the data you wish to plot
#' @param x.axis NO DEFAULT. X axis
#' @param y.axis NO DEFAULT. Y axis
#'
#' @param plot.by NO DEFAULT. A list or vector of column names (markers) to plot.
#' @param density.plot DEFAULT = FALSE. Logical. Creates a density plot (using 'ggpointdensity' package).
#' @param align.xy.by NO DEFAULT. Align X and Y to a dataset
#' @param align.col.by NO DEFAULT. Align colour to a dataset
#'
#' @param colours DEFAULTS to 'spectral'. What colour scheme do you want to use. Only used if type = 'colour', ignored if type = 'factor'. Can be 'jet', 'spectral', 'viridis', 'inferno', or 'magma'.
#' @param figure.title DEFAULTS to 'Marker expression'. Figure title
#' @param dot.size DEFAULTS to 1. Size of dots
#' @param col.min.threshold DEFAULTS to 0.01. Minimum threshold for colour scale.
#' @param col.max.threshold DEFAULTS to 0.995. Maximum threshold for colour scale.
#' @param path DEFAULTS TO getwd() -- i.e. the current working directory. Path to the desired output directory
#' @param plot.width DEFAULTS to 9.
#' @param plot.height DEFAULTS to 7.
#' @param blank.axis DEFAULT = FALSE. Logical, do you want a minimalist graph?
#' @param save.each.plot DEFAULT = FALSE. Do you want to save each plot?
#'
#' @usage make.multi.marker.plot(dat, x.axis, y.axis, plot.by, density.plot, align.xy.by, align.col.by, colours, figure.title, dot.size, col.min.threshold, col.max.threshold, path, plot.width,  plot.height, blank.axis, save.each.plot)
#'
#' @examples
#' # Create grid of plots on demonstration data
#' Spectre::make.multi.marker.plot(dat = Spectre::demo.umap,
#'                                 x.axis <- "UMAP_42_X",
#'                                 y.axis <- "UMAP_42_Y",
#'                                 plot.by <- c("AF700.CD45", "APCCy7.CD48", "BV605.Ly6C", "DL800.SA.Bio.Ly6G"),
#'                                 )
#'
#' @author 
#' Thomas Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#' @export

make.multi.marker.plot <- function(dat,
                               x.axis, # X axis for all plots
                               y.axis, # Y axis for all plots
                               plot.by, # Column -- where each unique value is plotted separately (i.e. plot by sample, etc)
                               density.plot = FALSE,
                               align.xy.by, # alignment for X and Y
                               align.col.by, # alignment for colours
                               colours = 'spectral',
                               figure.title = 'Marker expression',
                               dot.size = 1,
                               col.min.threshold = 0.01,
                               col.max.threshold = 0.995,
                               path = getwd(),
                               plot.width = 9,
                               plot.height = 7,
                               blank.axis = FALSE,
                               save.each.plot = FALSE
                              )
{

  ### Test data

      # library(Spectre)
      # library('plyr')
      # library('data.table')
      # library('tidyr') # for spread
      # library('rstudioapi')
      # library('ggplot2')
      # library('scales')
      # library('colorRamps')
      # library('ggthemes')
      # library('RColorBrewer')
      # library("gridExtra")
      #
      # library(data.table)
      #
      # ### COLOURS
      #
      # setwd("/Users/Tom/Desktop")
      # getwd()

      # dat <- Spectre::demo.umap
      # x.axis <- "UMAP_42_X"
      # y.axis <- "UMAP_42_Y"
      #
      # plot.by <- c("AF700.CD45", "APCCy7.CD48", "BV605.Ly6C", "DL800.SA.Bio.Ly6G")
      #
      # align.xy.by <- Spectre::demo.umap
      # align.col.by <- Spectre::demo.umap
      #
      # colours = 'spectral'
      # figure.title = 'Marker expression'
      # dot.size = 1
      # col.min.threshold = 0.01
      # col.max.threshold = 0.995
      # path = getwd()
      # plot.width = 9
      # plot.height = 7
      # blank.axis = FALSE
      # save.each.plot = FALSE

  ## Check that necessary packages are installed
  if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
  if(!is.element('ggplot2', installed.packages()[,1])) stop('ggplot2 is required but not installed')
  if(!is.element('scales', installed.packages()[,1])) stop('scales is required but not installed')
  if(!is.element('colorRamps', installed.packages()[,1])) stop('colorRamps is required but not installed')
  if(!is.element('ggthemes', installed.packages()[,1])) stop('ggthemes is required but not installed')
  if(!is.element('RColorBrewer', installed.packages()[,1])) stop('RColorBrewer is required but not installed')
  if(!is.element('ggpointdensity', installed.packages()[,1])) stop('ggpointdensity is required but not installed')

  ## Require packages
  require(Spectre)
  require(ggplot2)
  require(scales)
  require(colorRamps)
  require(ggthemes)
  require(RColorBrewer)
  require(ggpointdensity)

  ### Create each plot
  title <- figure.title

  plots <- list()
  to.plot <- plot.by

  ## Loop
    for(i in to.plot){
      plots[[i]] <- Spectre::make.colour.plot(dat = dat, #instead, use d[d[$Sample][plot.by] == to.plot[i], ]
                                x.axis = x.axis,
                                y.axis = y.axis,
                                col.axis = i,
                                title = i,
                                align.xy.by = dat,
                                align.col.by = dat,
                                colours = colours,
                                col.min.threshold = col.min.threshold,
                                col.max.threshold = col.max.threshold,
                                dot.size = dot.size,
                                path = path,
                                plot.width = plot.width,
                                plot.height = plot.height,
                                blank.axis = blank.axis,
                                save.to.disk = save.each.plot
                                )
    }

  ## Add density plot
  if(density.plot == TRUE){
    plots[[length(to.plot) + 1]] <- Spectre::make.density.plot(dat = dat,
                                                 x.axis = x.axis,
                                                 y.axis = y.axis,
                                                 colours = "viridis",
                                                 dot.size = dot.size,
                                                 align.xy.by = dat,
                                                 align.col.by = dat,
                                                 save.to.disk = save.each.plot,
                                                 path = path,
                                                 plot.width = plot.width,
                                                 plot.height = plot.height,
                                                 blank.axis = blank.axis
                                                 )
    # Rename density data set
    names(plots)[length(to.plot) + 1] <- "Density"
  }


  ### Arrange plots in grd

  if(length(plots) < 4){
    num.cols <- length(plots)
    num.rows <- 1
  }

  if(length(plots) == 4){
    num.cols <- 4
    num.rows <- 1
  }

  if(length(plots) > 4){
    num.cols <- 4
    num.rows <- length(plots)/4
    num.rows <- ceiling(num.rows)
  }

  gp <- gridExtra::grid.arrange(grobs = plots,
                     ncol = num.cols,
                     nrow = num.rows
                     #top = figure.title,
                     #top = textGrob(figure.title,gp=gpar(fontsize=20,font=3)),
  )
  #top = figure.title) #top = "Main Title" -- need to fix size issues

  ### Save to disk
  if (save.each.plot == TRUE) {
    # width = 9 per graph (4 graphs across max, so 4 cols max)
    wdth <- num.cols*plot.width
    
    # height = 7 per graph
    hght <- num.rows*plot.height
    
    ggplot2::ggsave(filename = paste0(figure.title, " - plotted on ", x.axis, " by ", y.axis, ".png"),
                    plot = gp,
                    path = path,
                    width = wdth,
                    height = hght,
                    limitsize = FALSE)
    
    ### Message
    
    if(exists(x = "gp")){
      message(paste0("Check your working directory for a new .png called ", "'",figure.title,".png'"))
    }
  }


  if(exists(x = "gp") == FALSE){
    message(paste0("Grid image was not created successfully"))
  }

}

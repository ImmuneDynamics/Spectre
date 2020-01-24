#' multi.marker.plot
#'
#' @param d NO DEFAULT. A data frame containing all the data you wish to plot
#' @param x.axis NO DEFAULT. X axis
#' @param y.axis NO DEFAULT. Y axis
#'
#' @param plot.by NO DEFAULT. A list or vector of column names (markers) to plot.
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
#'
#' This function allows you to create a grid of plots, where the cells are subsetted by a certain factor (e.g. one sample per plot). These can then be coloured by a marker or by another factor (e.g Group).
#'
#' @usage multi.marker.plot(x, ...)
#'
#' @export

multi.marker.plot <- function(d,
                               x.axis, # X axis for all plots
                               y.axis, # Y axis for all plots

                               plot.by, # Column -- where each unique value is plotted separately (i.e. plot by sample, etc)

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

  # d <- demo.umap
  # x.axis <- "UMAP_42_X"
  # y.axis <- "UMAP_42_Y"
  #
  # plot.by <- c("AF700.CD45", "APCCy7.CD48", "BV605.Ly6C", "DL800.SA.Bio.Ly6G")
  #
  # align.xy.by <- demo.umap
  # align.col.by <- demo.umap
  #
  # colours = 'spectral'
  # figure.title = 'Marker expression'
  # dot.size = 1
  # col.min.threshold = 0.01
  # col.max.threshold = 0.995
  # path = getwd()
  # plot.width = 9
  # plot.height = 7

  ### Create each plot
  plots <- list()
  to.plot <- plot.by

  ## Loop
    for(i in to.plot){
      plots[[i]] <- colour.plot(d = d, #instead, use d[d[$Sample][plot.by] == to.plot[i], ]
                                x.axis = x.axis,
                                y.axis = y.axis,
                                col.axis = i,
                                title = i,
                                align.xy.by = d,
                                align.col.by = d,
                                colours = colours,
                                col.min.threshold = col.min.threshold,
                                col.max.threshold = col.max.threshold,
                                dot.size = dot.size,
                                path = path,
                                plot.width = plot.width,
                                plot.height = plot.height,
                                save.to.disk = save.each.plot
                                )
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

  gp <- grid.arrange(grobs = plots,
                     ncol = num.cols,
                     nrow = num.rows,
                     #top = figure.title,
                     #top = textGrob(figure.title,gp=gpar(fontsize=20,font=3)),
  )
  #top = figure.title) #top = "Main Title" -- need to fix size issues

  ### Save to disk

  # width = 9 per graph (4 graphs across max, so 4 cols max)
  wdth <- num.cols*plot.width

  # height = 7 per graph
  hght <- num.rows*plot.height

  ggsave(filename = paste0(figure.title, ".png"),
         plot = gp,
         path = path,
         width = wdth,
         height = hght,
         limitsize = FALSE)

  ### Message

  if(exists(x = "gp")){
    print(paste0("Check your working directory for a new .png called ", "'",figure.title,".png'"))
  }

  if(exists(x = "gp") == FALSE){
    print(paste0("Grid image was not created successfully"))
  }

}

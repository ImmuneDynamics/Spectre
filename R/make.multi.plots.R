#' make.multi.plot
#'
#' @param dat NO DEFAULT. A data frame containing all the data you wish to plot
#' @param x.axis NO DEFAULT. X axis
#' @param y.axis NO DEFAULT. Y axis
#' @param col.axis NO DEFAULT. Colour axis
#' @param type NO DEFAULT. Specify whether col.axis is continious ('colour', for markers) or a factor ('factor', for clusters, sample, groups etc), 'labelled.factor' or 'density' (for density plots - note figure.title must be entered).
#'
#' @param plot.by NO DEFAULT. Name of the FACTORS that the plots will be divided by (e.g. Sample, Group, Batch etc).
#' @param align.xy.by NO DEFAULT. Align X and Y to a dataset
#' @param align.col.by NO DEFAULT. Align colour to a dataset
#'
#' @param colours DEFAULTS to 'spectral'. What colour scheme do you want to use. Only used if type = 'colour', ignored if type = 'factor'. Can be 'jet', 'spectral', 'viridis', 'inferno', or 'magma'.
#' @param figure.title DEFAULTS to paste0(plot.by, " ", col.axis). Figure title. NOTE, CURRENTLY DISABLED.
#' @param dot.size DEFAULTS to 1. Size of dots
#' @param col.min.threshold DEFAULTS to 0.01. Minimum threshold for colour scale.
#' @param col.max.threshold DEFAULTS to 0.995. Maximum threshold for colour scale.
#' @param path DEFAULTS TO getwd() -- i.e. the current working directory. Path to the desired output directory
#' @param plot.width DEFAULTS to 9.
#' @param plot.height DEFAULTS to 7.
#' @param blank.axis DEFAULT = FALSE. Logical, do you want a minimalist graph?
#' @param save.each.plot DEFAULT = FALSE. Do you want to save each plot?
#'
#' This function allows you to create a grid of plots, where the cells are subsetted by a certain factor (e.g. one sample per plot). These can then be coloured by a marker or by another factor (e.g Group).
#'
#' @usage make.multi.plot(dat, x.axis, y.axis, col.axis, type, plot.by, align.xy.by, align.col.by, colours, figure.title, dot.size, col.min.threshold, col.max.threshold, path, plot.width,  plot.height, blank.axis, save.each.plot, ...)
#'
#' @export

make.multi.plot <- function(dat,
                       x.axis, # X axis for all plots
                       y.axis, # Y axis for all plots
                       col.axis,
                       type, # colour, factor, labelled.factor

                       plot.by, # Column -- where each unique value is plotted separately (i.e. plot by sample, etc)

                       align.xy.by, # alignment for X and Y
                       align.col.by, # alignment for colours

                       colours = 'spectral',
                       figure.title = paste0(col.axis, " by ", plot.by),
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
      # ### FACTORS
      #
      # setwd("/Users/Tom/Desktop")
      # getwd()
      #
      # dat <- Spectre::demo.umap
      # x.axis = "UMAP_42_X"
      # y.axis = "UMAP_42_Y"
      # col.axis = "Group"
      # type = "factor"
      #
      # plot.by = "Sample"
      # align.xy.by = dat
      # align.col.by = dat
      #
      # colour = "magma"
      # figure.title = paste0("By ", plot.by, " - ", col.axis)
      # dot.size = 3
      # path = getwd()
      # plot.width = 9
      # plot.height = 7
      #
      # ### COLOURS
      #
      # setwd("/Users/Tom/Desktop")
      # getwd()
      #
      # dat <- Spectre::demo.umap
      # x.axis = "UMAP_42_X"
      # y.axis = "UMAP_42_Y"
      # col.axis = "AF700.CD45"
      # type = "colour"
      #
      # plot.by = "Sample"
      # align.xy.by = dat
      # align.col.by = dat
      #
      # colours = "magma"
      # figure.title = paste0("By ", plot.by, " - ", col.axis)
      # dot.size = 3
      # col.min.threshold = 0.01
      # col.max.threshold = 0.995
      # path = getwd()
      # plot.width = 9
      # plot.height = 7

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
      to.plot <- list()
      plots <- list()

      title <- figure.title

      to.plot <- unique(dat[[plot.by]])

      ## For type = 'factor'
      if(type == 'factor'){
        for(i in c(1:length(to.plot))){
          to.plot[i]
          plots[[i]] <- Spectre::make.factor.plot(dat = subset(dat, dat[[plot.by]] == to.plot[i]), #instead, use d[d[$Sample][plot.by] == to.plot[i], ]
                                    x.axis = x.axis,
                                    y.axis = y.axis,
                                    col.axis = col.axis,
                                    title = to.plot[i],
                                    align.xy.by = dat,
                                    align.col.by = dat,
                                    dot.size = dot.size,
                                    path = path,
                                    plot.width = plot.width,
                                    plot.height = plot.height,
                                    blank.axis = blank.axis,
                                    save.to.disk = save.each.plot)
        }
      }

      ## For type = 'factor'
      if(type == 'labelled.factor'){
        for(i in c(1:length(to.plot))){
          to.plot[i]
          plots[[i]] <- Spectre::make.factor.plot(dat = subset(dat, dat[[plot.by]] == to.plot[i]), #instead, use d[d[$Sample][plot.by] == to.plot[i], ]
                                                  x.axis = x.axis,
                                                  y.axis = y.axis,
                                                  col.axis = col.axis,
                                                  title = to.plot[i],
                                                  align.xy.by = dat,
                                                  align.col.by = dat,
                                                  dot.size = dot.size,
                                                  path = path,
                                                  plot.width = plot.width,
                                                  plot.height = plot.height,
                                                  blank.axis = blank.axis,
                                                  save.to.disk = save.each.plot, add.label = TRUE)
        }
      }


      ## For type = 'colour'
      if(type == 'colour'){
          for(i in c(1:length(to.plot))){
            to.plot[i]
            plots[[i]] <- Spectre::make.colour.plot(dat = subset(dat, dat[[plot.by]] == to.plot[i]), #instead, use d[d[$Sample][plot.by] == to.plot[i], ]
                                      x.axis = x.axis,
                                      y.axis = y.axis,
                                      col.axis = col.axis,
                                      title = to.plot[i],
                                      align.xy.by = dat,
                                      align.col.by = dat,
                                      colours = colours,
                                      col.min.threshold = col.min.threshold,
                                      col.max.threshold = col.max.threshold,
                                      dot.size = dot.size,
                                      blank.axis = blank.axis)
          }
      }

      if(type == 'density'){
          for (i in c(1:length(to.plot))){
            to.plot[i]
            plots[[i]] <- Spectre::make.density.plot(dat = subset(dat, dat[[plot.by]] == to.plot[i]), #instead, use d[d[$Sample][plot.by] == to.plot[i], ]
                                                     x.axis = x.axis,
                                                     y.axis = y.axis,
                                                     title = to.plot[i],
                                                     colours = "viridis",
                                                     dot.size = dot.size,
                                                     align.xy.by = dat,
                                                     align.col.by = dat,
                                                     save.to.disk = save.each.plot,
                                                     path = path,
                                                     plot.width = plot.width,
                                                     plot.height = plot.height,
                                                     blank.axis = blank.axis)
          }
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

          # width = 9 per graph (4 graphs across max, so 4 cols max)
          wdth <- num.cols*plot.width

          # height = 7 per graph
          hght <- num.rows*plot.height

          ggplot2::ggsave(filename = paste0("Multi plot - ", figure.title, " - plotted on ", x.axis, " by ", y.axis, ".png"),
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

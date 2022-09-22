#' Make multiple plots for multiple columns and/or multiple sample/group/clusters
#'
#' Method to create multiple plots for each marker.
#' This function allows you to create a grid of plots, where the cells are plotted by a series of columns, and/or subsetted by a certain factor (e.g. one sample per plot).
#' Makes use of Spectre functions make.colour.plot and make.density.plot.
#'
#' @param dat NO DEFAULT. A data frame containing all the data you wish to plot
#' @param x.axis NO DEFAULT. X axis
#' @param y.axis NO DEFAULT. Y axis
#' @param plot.by NO DEFAULT. A vector of character names for the columns you wish to plot.
#'
#' @param divide.by DEFAULT = NULL. Here you can specify a character name of a column you wish to use to divide up the dataset.
#' @param add.density DEFAULT = FALSE. Can specify to add a density plot at the end the series of colour plots
#'
#' @param col.type DEFAULT = "continuous". Can also be "factor".
#' @param hex DEFAULT = FALSE. Whether to split the data into bins and show the average expression of the bin. Currently does not work with density plots, only for those features in the plot.by.
#' @param hex.bins DEFAULT = 30. Number of bins to split into. Only used if hex is TRUE.
#' @param figure.title DEFAULT = "Multi plot". Also used as the prefix for the saved file name.
#' @param global.xy DEFAULT = TRUE. Defines the limits for the X and Y based on the whole dataset. If FALSE, then each plot X & Y limits scale individually.
#' @param global.col DEFAULT = TRUE. Defines the limits for the colour axis based on the whole dataset. If FALSE, then each plot colour limit scales individually.
#' @param colours DEFAULTS to 'spectral'. What colour scheme do you want to use. Only used if type = 'colour', ignored if type = 'factor'. Can be 'jet', 'spectral', 'viridis', 'inferno', 'magma', or "BuPu".
#' @param dot.size DEFAULT = 1. Numeric. Size of the dots.
#' @param col.min.threshold DEFAULT = 0.01. Numeric. Define minimum threshold for colour scale. Values below this limit will be coloured as the chosen minimum threshold.
#' @param col.max.threshold DEFAULT = 0.995 Numeric. Define maximum threshold for colour scale. Values above this limit will be coloured as the chosen maximum threshold.
#' @param path DEFAULT = getwd() -- i.e. the current working directory. Path to the desired output directory
#' @param plot.width DEFAULT = 9.
#' @param plot.height DEFAULT = 7.
#' @param blank.axis DEFAULT = FALSE. Logical. Do you want a minimalist graph?
#' @param save.each.plot DEFAULT = FALSE. Logical. Do you want to save each plot?
#'
#' @usage make.multi.plot(dat, x.axis, y.axis, plot.by, divide.by, add.density, col.type, figure.title, align.xy.by, align.col.by, colours, dot.size, col.min.threshold, col.max.threshold, path, plot.width, plot.height, blank.axis, save.each.plot)
#'
#' @examples
#' # Create grid of plots on demonstration data
#' Spectre::make.multi.plot(dat = as.data.table(Spectre::demo.umap),
#'                           x.axis = "UMAP_42_X",
#'                           y.axis = "UMAP_42_Y",
#'                           plot.by = c("BV605.Ly6C", "BUV737.B220", "AF700.CD45"))
#'
#' @author
#' Thomas Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#'
#' @export

# align.xy.by DEFAULT = dat. Align X and Y to a dataset. By default it will be based on the total dataset.
# align.col.by DEFAULT = dat. Align colour to a dataset. By default it will be based on the total dataset.

make.multi.plot <- function(dat,
                            x.axis,
                            y.axis,
                            plot.by, # vector of column names -- one colour plot will be created for each
                            
                            divide.by = NULL,
                            add.density = FALSE,
                            
                            hex = FALSE,
                            hex.bins = 30,
                            
                            col.type = "continuous",
                            figure.title = 'Multi plot',
                            
                            global.xy = TRUE,
                            global.col = TRUE,
                            
                            align.xy.by = dat, # alignment for X and Y
                            align.col.by = dat, # alignment for colours
                            
                            colours = 'spectral',
                            dot.size = 1,
                            col.min.threshold = 0.01,
                            col.max.threshold = 0.995,
                            path = getwd(),
                            plot.width = 9, # each plot
                            plot.height = 7, # each plot
                            blank.axis = FALSE,
                            save.each.plot = FALSE){
  
  ### Check packages
  if(!is.element('ggplot2', installed.packages()[,1])) stop('ggplot2 is required but not installed')
  if(!is.element('scales', installed.packages()[,1])) stop('scales is required but not installed')
  if(!is.element('colorRamps', installed.packages()[,1])) stop('colorRamps is required but not installed')
  if(!is.element('ggthemes', installed.packages()[,1])) stop('ggthemes is required but not installed')
  if(!is.element('RColorBrewer', installed.packages()[,1])) stop('RColorBrewer is required but not installed')
  if(!is.element('ggpointdensity', installed.packages()[,1])) stop('ggpointdensity is required but not installed')
  
  ### Load packages
  require(ggplot2)
  require(scales)
  require(colorRamps)
  require(ggthemes)
  require(RColorBrewer)
  require(ggpointdensity)
  
  ### Overall plot settings
  #title <- figure.title
  
  plots <- list()
  to.plot <- plot.by
  
  ### Divide by, if requested
  
  dat.list <- list()
  
  if(is.null(divide.by)){ # the complete dataset
    dat.list[["Only"]] <- dat
  }
  
  if(!is.null(divide.by)){ # divide dataset by sample/group etc
    to.divide <- unique(dat[[divide.by]])
    
    for(i in to.divide){
      rws <- dat[[divide.by]] == i
      dat.list[[i]] <- dat[rws,]
    }
  }
  
  
  ### Loop - plot each colour
  
  for(a in names(dat.list)){
    # a <- "Only"
    # a <- "Mock"
    
    plot.dat <- dat.list[[a]]
    
    ## Create plots
    for(i in to.plot){
      
      if(length(dat.list) == 1){
        plot.nme <- paste0(i)
      }
      
      if(length(dat.list) > 1){
        plot.nme <- paste0(a, " - ", i)
      }
      
      if(global.xy == FALSE){
        align.xy.by <- plot.dat
      }
      
      if(global.col == FALSE){
        align.col.by <- plot.dat
      }
      
      plots[[plot.nme]] <- make.colour.plot(dat = plot.dat, #instead, use d[d[$Sample][plot.by] == to.plot[i], ] ## Spectre::
                                            x.axis = x.axis,
                                            y.axis = y.axis,
                                            col.axis = i,
                                            col.type = col.type,
                                            hex = hex,
                                            hex.bins = hex.bins,
                                            title = plot.nme,
                                            align.xy.by = align.xy.by, ###
                                            align.col.by = align.col.by, ###
                                            colours = colours,
                                            col.min.threshold = col.min.threshold,
                                            col.max.threshold = col.max.threshold,
                                            dot.size = dot.size,
                                            path = path,
                                            plot.width = plot.width,
                                            plot.height = plot.height,
                                            blank.axis = blank.axis,
                                            save.to.disk = save.each.plot)
    }
    
    ## Add density plot (if desired)
    if(length(dat.list) == 1){
      if(add.density == TRUE){
        
        
        if(global.xy == TRUE){
          align.xy.by <- dat
        } else {
          align.xy.by <- plot.dat
        }
        
        plots[[length(plots) + 1]] <- make.colour.plot(dat = plot.dat,
                                                       x.axis = x.axis,
                                                       y.axis = y.axis,
                                                       colours = colours,
                                                       dot.size = dot.size,
                                                       title = "Density",
                                                       align.xy.by = align.xy.by,
                                                       save.to.disk = save.each.plot,
                                                       path = path,
                                                       plot.width = plot.width,
                                                       plot.height = plot.height,
                                                       blank.axis = blank.axis)
        names(plots)[length(plots)] <- paste0("Density")
      }
    }
    
    if(length(dat.list) > 1){
      if(add.density == TRUE){
        
        if(global.xy == TRUE){
          align.xy.by <- dat
        } else {
          align.xy.by <- plot.dat
        }
        
        plots[[length(plots) + 1]] <- make.colour.plot(dat = plot.dat,
                                                       x.axis = x.axis,
                                                       y.axis = y.axis,
                                                       colours = colours,
                                                       dot.size = dot.size,
                                                       title = paste0(a, "_Density"),
                                                       align.xy.by = align.xy.by,
                                                       save.to.disk = save.each.plot,
                                                       path = path,
                                                       plot.width = plot.width,
                                                       plot.height = plot.height,
                                                       blank.axis = blank.axis)
        names(plots)[length(plots)] <- paste0(a, "_Density")
      }
    }
  }
  
  
  ### Arrange plots in grid
  
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
                                nrow = num.rows)
  #top = figure.title,
  #top = textGrob(figure.title,gp=gpar(fontsize=20,font=3)),
  #top = figure.title) #top = "Main Title" -- need to fix size issues
  
  ### Save to disk
  # if (save.each.plot == TRUE) {
  #   message("The option to save each plot individually is not currently available")
  # }
  
  # width = 9 per graph (4 graphs across max, so 4 cols max)
  wdth <- num.cols*plot.width
  
  # height = 7 per graph
  hght <- num.rows*plot.height
  
  if(add.density == TRUE){
    figure.title <- paste0(figure.title, " plus density")
  }
  
  if(!is.null(divide.by)){
    figure.title <- paste0(figure.title, " divided by ", divide.by)
  }
  
  ggplot2::ggsave(filename = paste0(figure.title, " - plotted on ", x.axis, " by ", y.axis, ".png"),
                  plot = gp,
                  path = path,
                  width = wdth,
                  height = hght,
                  limitsize = FALSE)
  
  if(exists(x = "gp")){
    message(paste0("Check your working directory for a new .png called ", "'",figure.title,".png'"))
  }
  
  ### Message if grid was NOT created
  
  if(exists(x = "gp") == FALSE){
    message(paste0("Grid image was not created successfully"))
  }
}

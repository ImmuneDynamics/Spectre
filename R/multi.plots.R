#' multi.plot
#'
#' @param d NO DEFAULT. A data frame containing all the data you wish to plot
#' @param type DEFAULTS to 'factor' to create factor plots. Will add 'colour' option to create colour plots in next update.
#' @param x.axis X axis 
#' @param y.axis Y axis
#' @param plot.by What factor should the plots be divided by
#' @param align.xy.by Align X and Y to a dataset
#' @param align.col.by Align colour to a dataset
#' @param colour What feature should the plot be coloured by
#' @param figure.title Figure title
#' @param dot.size Size of dots
#'
#' @usage multi.plot(x, ...)
#'
#' @export

multi.plot <- function(d, 
                       type = 'factor', # colour, factor
                       x.axis, # X axis for all plots
                       y.axis, # Y axis for all plots
                       plot.by, # Column -- where each unique value is plotted separately (i.e. plot by sample, etc)
                       align.xy.by, # alignment for X and Y 
                       align.col.by, # alignment for colours
                       colour,
                       figure.title,
                       dot.size
                       )
{
  
  ### Test data
  
      setwd("/Users/Tom/Desktop")
      getwd()
  
      d <- demo.umap
      type = "factor"
      x.axis = "UMAP_42_X"
      y.axis = "UMAP_42_Y" 
      plot.by = "Sample"
      align.xy.by = d
      align.col.by = d
      colour = "magma"
      figure.title = "By sample"
      dot.size = 1
  
  ### Create each plot
  
      ## For type = 'factor'
          to.plot <- list()
          plots <- list()
          
          to.plot <- unique(d[[plot.by]])
          
          for(i in c(1:length(to.plot))){
            to.plot[i]
            plots[[i]] <- factor.plot(d = subset(d, d[[plot.by]] == to.plot[i]), #instead, use d[d[$Sample][plot.by] == to.plot[i], ]
                             x.axis = x.axis,
                             y.axis = y.axis, 
                             col.axis = plot.by,
                             title = to.plot[i],
                             align.xy.by = d,
                             align.col.by = d,
                             dot.size = dot.size)
            }
      
      ## For type = 'colour'
          # plots <- list()
          # plots <- unique(d[[plot.by]])
          # 
          # for(i in plots){
          #   p <- colour.plot(d = d,
          #                   x.axis = x.axis,
          #                   y.axis = y.axis, 
          #                   col.axis = plot.by,
          #                   colour = colour,
          #                   title = plot.by)
          #   
          #   
          # }
  
  
  
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
          }
          
          gp <- grid.arrange(grobs = plots, 
                             ncol = num.cols, 
                             nrow = num.rows)
                             #top = figure.title) #top = "Main Title" -- need to fix size issues
  
  ### Save to disk
          
          # width = 9 per graph (4 graphs across max, so 4 cols max)
          wdth <- num.cols*9
          
          # height = 7 per graph
          hght <- num.rows*7
          
          ggsave(filename = paste0(figure.title, ".png"), 
                 plot = gp, 
                 path = getwd(), 
                 width = wdth, 
                 height = hght)
          
  ### Message
          
          if(exists(x = "gp")){
            print(paste0("Check your working directory for a new .png called ", "'",figure.title,".png'"))
          }
          
          if(exists(x = "gp") == FALSE){
            print(paste0("Grid image was not created successfully"))
          }
      
}

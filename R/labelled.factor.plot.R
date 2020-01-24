#' labelled.factor.plot - Create a factor plot
#'
#' @usage labelled.factor.plot(x, ...)
#'
#' @param d data.frame. Input sample. No default.
#' @param a.axis Character. Column for X axis. No default.
#' @param y.axis Character. Column for Y axis. No default.
#' @param col.axis Character. Column for colour. No default.
#' @param title Character. Title for the plot. Defaults to the name of the colour axis.
#' @param dot.size Numeric. Size of the dots. Defaults to 1.
#' @param align.xy.by data.frame. Sample to use to determine minimum and maximum X and Y axis values. No default.
#' @param align.col.by data.frame. Sample to use to determine minimum and maximum colour values. No default.
#' @param nudge_x DEFAULTS to 0.5. X-nudge distance for labels.
#' @param nudge_y DEFAULTS to 0.5. X-nudge distance for labels.
#'
#' This function allows you to plot cells on an XY plot, coloured by some factor (group, cluster, etc) with labels included.
#' factor.plot is a wrapper around the heatmap.2 function from gplots.
#'
#' @export

labelled.factor.plot <- function(d,
                                 x.axis, # "UMAP1"
                                 y.axis, # "UMAP2"
                                 col.axis, # "BV605.Ly6C"
                                 title = col.axis,
                                 dot.size = 1,
                                 align.xy.by = NULL,
                                 align.col.by = NULL,
                                 nudge_x = 0.5,
                                 nudge_y = 0.5,
                                 save.to.disk = TRUE,
                                 path = getwd(),
                                 plot.width = 9,
                                 plot.height = 7){
  ## TESTING
      # install.packages('ggrepel')
      # library('ggrepel')
      #
      # d = demo.umap
      # x.axis = "UMAP_42_X"
      # y.axis = "UMAP_42_Y"
      # col.axis = "FlowSOM_metacluster" # col.axis = "Sample"
      # title = "Cluster"
      # dot.size = 0.5
      # align.xy.by = demo.umap
      # align.col.by = demo.umap


  ## X, Y, col lims

        # X AXIS
        if(is.null(align.xy.by) == TRUE){
          Xmax <- max(d[[x.axis]])
          Xmin <- min(d[[x.axis]])
        }

        if(is.null(align.xy.by) == FALSE){
          Xmax <- max(align.xy.by[[x.axis]])
          Xmin <- min(align.xy.by[[x.axis]])
        }


        # Y AXIS
        if(is.null(align.xy.by) == TRUE){
          Ymax <- max(d[[y.axis]])
          Ymin <- min(d[[y.axis]])
        }

        if(is.null(align.xy.by) == FALSE){
          Ymax <- max(align.xy.by[[y.axis]])
          Ymin <- min(align.xy.by[[y.axis]])
        }

        # COLOUR
        if(is.null(align.col.by) == TRUE){
          #ColrMin <- min(d[[col.axis]])
          #ColrMax <- max(d[[col.axis]])

          colRange <- unique(d[[col.axis]])
          colRange <- colRange[order(colRange)]
          colRange <- as.character(colRange)
        }

        if(is.null(align.col.by) == FALSE){
          #ColrMin <- min(align.col.by[[col.axis]])
          #ColrMax <- max(align.col.by[[col.axis]])

          colRange <- unique(align.col.by[[col.axis]])
          colRange <- colRange[order(colRange)]
          colRange <- as.character(colRange)
        }


  ## Prep centroid coordinates


    if(is.numeric(d[[col.axis]])){
      centroidsDf <- data.frame(
        centroidX = tapply(d[[x.axis]], d[[col.axis]], median), # median
        centroidY = tapply(d[[y.axis]], d[[col.axis]], median),
        centroidCol = tapply(d[[col.axis]], d[[col.axis]], median))
    }

    if(!is.numeric(d[[col.axis]])){
      labels <- sort(unique(d[[col.axis]]))

      centroidsDf <- data.frame(
        centroidX = tapply(d[[x.axis]], d[[col.axis]], median), # median
        centroidY = tapply(d[[y.axis]], d[[col.axis]], median),
        centroidCol = labels)
    }

    print(ggplot(data = d, aes(x = d[[x.axis]], y = d[[y.axis]], colour = as.factor(d[[col.axis]]))) +
            geom_point(size = dot.size)+ # 2 for large # 0.5 for small
            #scale_colour_gradientn(colours = colour.scheme(50)) +
            #scale_colour_manual(name = "FileName", values = c(colour.scheme(length(filenames)))) +
            ggtitle(title) +
            #labs(colour = col.axis)+
            xlab(x.axis)+
            ylab(y.axis)+
            xlim(Xmin, Xmax) +
            ylim(Ymin, Ymax) +
            lims(colour = colRange) +

            theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # change 'colour' to black for informative axis
                  #axis.line=element_blank(),
                  #axis.text.x=element_blank(),
                  #axis.text.y=element_blank(),
                  #axis.ticks=element_blank(),
                  axis.title.x=element_text(color="Black", face="bold", size=18),
                  axis.title.y=element_text(color="Black", face="bold", size=18),
                  #panel.grid.major = element_blank(),
                  #panel.background=element_blank(),
                  #panel.border=element_blank(),
                  #panel.grid.minor=element_blank(),
                  #plot.background=element_blank(),
                  #legend.position = "right",
                  #legend.text=element_text(size=12), # large = 30 # small = 8
                  #legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
                  #legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
                  legend.title=element_blank(),
                  plot.title = element_text(color="Black", face="bold", size=22, hjust=0) # size 70 for large, # 18 for small
            )+

            geom_point(data = centroidsDf,                     # separate data.frame
                       aes(x = centroidX, y = centroidY),
                       col = "black",                          # notice "col" and "shape" are
                       #shape = 1,
                       size = 2) +

            geom_label(data = centroidsDf,
                       hjust = 0,
                       nudge_x = nudge_x,
                       nudge_y = nudge_y,
                       aes(x = centroidX, y = centroidY, label=centroidCol, alpha = 0.5),
                       col = "black",
                       fontface = "bold") +

            guides(alpha = "none")
    )

    if(save.to.disk == TRUE){
      ggsave(filename = paste0(title, ".png"),
             plot = p,
             path = path,
             width = plot.width,
             height = plot.height,
             limitsize = FALSE)
    }

}

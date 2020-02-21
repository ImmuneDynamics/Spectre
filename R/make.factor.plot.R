#' make.factor.plot - Create a dot plot (X vs Y) coloured by a selected factor column (e.g. clusters, samples)
#'
#' This function allows you to create a coloured XY plot where each cell is coloured by a selected column. Typically used to plot cells on tSNE1/2 or UMAP1/2 coloured by select factor, such as clusters, samples, groups, batches etc.
#'
#' @usage make.factor.plot(dat, x.axis, y.axis, col.axis)
#'
#' @param dat NO DEFAULT. data.table Input sample.
#' @param x.axis NO DEFAULT. Character. Column for X axis.
#' @param y.axis NO DEFAULT. Character. Column for Y axis.
#' @param col.axis NO DEFAULT. Character. Column for colour, should be a factor (e.g. cluster, group, sample, etc).
#' @param add.label DEFAULT = FALSE. Add's a label for the entries in 'col.axis' onto the plot -- useful for labelling clusters (e.g. 1, 2, 3 etc) or populations (T cell, B cell, etc).
#' @param nudge_x DEFAULT = 0.5. Distance to nudge the label from the point (on the x axis) when add.label = TRUE.
#' @param nudge_y DEFAULT = 0.5. Distance to nudge the label from the point (on the y axis) when add.label = TRUE.
#' @param title DEFAULT = col.axis. Character. Title for the plot.
#' @param col.min.threshold DEFAULT = 0.01. Numeric. Define minimum threshold for colour scale. Values below this limit will be coloured as the chosen minimum threshold.
#' @param col.max.threshold DEFAULT = 0.995 Numeric. Define maximum threshold for colour scale. Values above this limit will be coloured as the chosen maximum threshold.
#' @param dot.size DEFAULT = 1. Numeric. Size of the dots.
#' @param align.xy.by DEFAULT = dat. data.table Sample to use to determine minimum and maximum X and Y axis values.
#' @param align.col.by DEFAULT = dat. data.table. Sample to use to determine minimum and maximum colour values.
#' @param blank.axis DEFAULT = FALSE. Logical, do you want a minimalist graph?
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
#' make.factor.plot(dat = demo.umap, x.axis = "UMAP_X", y.axis = "UMAP_Y", col.axis = "FlowSOM_metacluster")
#'
#' @export

make.factor.plot <- function(dat,
                              x.axis,
                              y.axis,
                              col.axis,
                              add.label = FALSE,
                              nudge_x = 0.5,
                              nudge_y = 0.5,
                              title = col.axis,
                              dot.size = 1,
                              align.xy.by = dat,
                              align.col.by = dat,
                              save.to.disk = TRUE,
                              path = getwd(),
                              plot.width = 9,
                              plot.height = 7,
                              blank.axis = FALSE)
{

    ## TESTING
        # dat = demo.umap
        # x.axis = "UMAP_42_X"
        # y.axis = "UMAP_42_Y"
        # col.axis = "FlowSOM_metacluster"
        # add.label = TRUE
        # nudge_x = 0.5
        # nudge_y = 0.5
        # title = col.axis
        # dot.size = 1
        # align.xy.by = dat
        # align.col.by = dat
        # save.to.disk = TRUE
        # path = getwd()
        # plot.width = 9
        # plot.height = 7
        # blank.axis = FALSE

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
              #ColrMin <- min(d[[col.axis]])
              #ColrMax <- max(d[[col.axis]])

              colRange <- unique(dat[[col.axis]])
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

    ## Generate and show ggplot

        p <- ggplot(data = dat, aes(x = dat[[x.axis]], y = dat[[y.axis]], colour = as.factor(dat[[col.axis]]))) +
                    geom_point(size = dot.size) +
                    ggtitle(title) +
                    xlab(x.axis) +
                    ylab(y.axis) +
                    xlim(Xmin, Xmax) +
                    ylim(Ymin, Ymax) +
                    lims(colour = colRange)

    ## Add some themes

        p <- p + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5),
                        axis.title.x=element_text(color="Black", face="bold", size=18),
                        axis.title.y=element_text(color="Black", face="bold", size=18),
                        legend.title=element_blank(),
                        plot.title = element_text(color="Black", face="bold", size=22, hjust=0) # size 70 for large, # 18 for small
                      )

    ## Add labels if desired

        if(add.label == TRUE){

          # Prepare centroids
              if(is.numeric(dat[[col.axis]])){

                centroidX = tapply(dat[[x.axis]], dat[[col.axis]], median) # median
                centroidY = tapply(dat[[y.axis]], dat[[col.axis]], median)
                centroidCol = tapply(dat[[col.axis]], dat[[col.axis]], median)

                centroidsDf <- data.frame(centroidX, centroidY, centroidCol)
              }

              if(!is.numeric(dat[[col.axis]])){
                labels <- sort(unique(dat[[col.axis]]))

                centroidsDf <- data.frame(
                  centroidX = tapply(dat[[x.axis]], dat[[col.axis]], median), # median
                  centroidY = tapply(dat[[y.axis]], dat[[col.axis]], median),
                  centroidCol = labels)
              }

          # Add labels

              p <- p + geom_point(data = centroidsDf,                     # separate data.frame
                                  aes(x = centroidsDf$centroidX, y = centroidsDf$centroidY),
                                  col = "black",                          # notice "col" and "shape" are
                                  #shape = 1,
                                  size = 2)

              p <- p + geom_label(data = centroidsDf,
                                  hjust = 0,
                                  nudge_x = nudge_x,
                                  nudge_y = nudge_y,
                                  aes(x = centroidsDf$centroidX, y = centroidsDf$centroidY, label=centroidsDf$centroidCol, alpha = 0.5),
                                  col = "black",
                                  fontface = "bold")

              p <- p + guides(alpha = "none")
        }


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

    ## Save to disk if desired

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

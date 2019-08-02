#################################################################################
###### CytoTools - heatmap options
#################################################################################

### 1.1 - Install packages (if not already installed)
if(!require('flowCore')) {install.packages('flowCore')}
if(!require('plyr')) {install.packages('plyr')}
if(!require('data.table')) {install.packages('data.table')}
if(!require('rstudioapi')) {install.packages('rstudioapi')}
if(!require('R.utils')) {install.packages('R.utils')}
if(!require('ggplot2')) {install.packages('ggplot2')}
if(!require('gplots')) {install.packages('gplots')}
if(!require('RColorBrewer')) {install.packages('RColorBrewer')}
if(!require('tidyr')) {install.packages('tidyr')}
if(!require('Biobase')) {install.packages('Biobase')}

### 1.2 Load packages
library('flowCore')
library('plyr')
library('data.table')
library('rstudioapi')
library('R.utils')
library('ggplot2')
library('gplots')
library('RColorBrewer')
library('tidyr')
library('Biobase')

### CytoTools::make.heatmap

    #' Function for generating heatmaps
    #'
    #' This function allows...
    #' @param file.loc What is the location of your files? ?? Defaults to TRUE ??
    #' @keywords cats
    #' @export
    #' @examples
    #' cat_function()

    make.heatmap <- function(x,
                             sample.col,
                             col.rmv,
                             plot.title,
                             do.transpose,
                             do.normalise,
                             do.dendrogram,
                             dendro.set,
                             row.clustering
                             n.row.groups,
                             col.clustering,
                             n.col.groups){

      ### 2.1 - Embed cluster or population name as row name
      heatmap.data <- x
      rownames(heatmap.data) <- t(heatmap.input[samp.col])
      heatmap.data

      heatmap.data <- heatmap.data[-c(col.rmv)] # remove columns
      heatmap.data

      ### 2.2 - Transpose (ONLY IF REQUIRED) -- the longest set (clusters or parameters) on x-axis -- by default MARKERS are columns, CLUSTERS are rows -- transpose to flip these defaults
      if(do.transpose == 1){
        heatmap.data.t <- as.data.frame(t(heatmap.data))
        heatmap.data <- heatmap.data.t
      }

      ### 2.3 - NORMALISE BY COLUMN (i.e. each column/parameter has a max of 1 and a minimum of 0) # This is optional, but allows for better comparison between markers
      if(do.normalise == 1){
          row.nam <- row.names(heatmap.data)

          norm.fun <- function(x) {(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) -min(x, na.rm=TRUE))}
          heatmap.data.norm <- as.data.frame(lapply(heatmap.data, norm.fun)) # by default, removes the names of each row
          max(heatmap.data.norm)
          max(heatmap.data.norm[,5])
          heatmap.data.norm <- as.matrix(heatmap.data.norm)
          heatmap.data <- heatmap.data.norm
          rownames(heatmap.data) <- row.nam # add row names back
      }

      # convert to matrix
      heatmap.data <- as.matrix(heatmap.data)


      ### 2.4 - Set up clustering

      if(do.dendrogram == 1){

        # set the custom distance and clustering functions, per your example
        hclustfunc <- function(x) hclust(x, method="complete")
        distfunc <- function(x) dist(x, method="euclidean")

        # perform clustering on rows and columns
        cl.row <- hclustfunc(distfunc(heatmap.data))
        cl.col <- hclustfunc(distfunc(t(heatmap.data)))

        # work out no cols and rows
        nrow(heatmap.data)
        ncol(heatmap.data)

        # extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
        gr.row <- cutree(cl.row, n.row.groups)
        gr.col <- cutree(cl.col, n.col.groups)

        col1 <- group.col(n.row.groups)
        col2 <- group.col(n.col.groups)
      }


      ### 2.6 - Set up fold-change or normal preferences

      if(is.foldchange == 1){
        map.colour <- fold.palette(31)
        sym.key <- FALSE # TRUE need if NOT creating own breaks
        sym.breaks <- TRUE
        my.breaks <- seq(min.range, max.range, length.out = 32)
      }

      if(is.foldchange == 0){
        map.colour <- colour.palette
        sym.key <- FALSE
        sym.breaks <- FALSE
        my.breaks <- seq(min(heatmap.data), max(heatmap.data), length.out = 32)
      }

      scale.set <- "none" # "none", "both", "row", "column"


      ### 2.7 - Plot heatmap and dendrograms

      dev.off()
      #tiff(file=paste0("Output_MFI/heatmap", ".tiff"), width = 2300, height = 1600, compression = "none") # , width = 9, height = 8)
      par(cex.main=1.3)
      par(xpd = TRUE)

      HeatMap <- heatmap.2(as.matrix(heatmap.data),  #### http://stanstrup.github.io/heatmaps.html # https://stackoverflow.com/questions/22278508/how-to-add-colsidecolors-on-heatmap-2-after-performing-bi-clustering-row-and-co
                           main=plot.title,
                           notecol= "black", # font colour of all cells to black
                           key=TRUE,
                           #keysize=0.75,#key.title = NA,
                           #key.par = list(cex=1),

                           hclust=hclustfunc, distfun=distfunc,

                           dendrogram =dendro.set, # "both", "row", "column" etc
                           Rowv = row.clustering, # turns on row clustering # = as.dendrogram(cluster)
                           Colv = col.clustering, # turns on column clustering # = reorder(as.dendrogram(cluster.row), 10:1) to flip around, or maybe dictage
                           RowSideColors = col1[gr.row],
                           ColSideColors = col2[gr.col], # as.character(clusters)

                           revC = FALSE, # default FALSE
                           symm=FALSE, # default FALSE
                           symkey= sym.key, # default FALSE, TRUE FOR SYMMETRY IN FOLD CHANGE SITUATIONS
                           symbreaks=sym.breaks, # default FALSE

                           breaks=my.breaks,

                           scale=scale.set, # "none" or "row" or "column". default is to scale by row, I think? # scale argument only for heat data, not dendrogram
                           trace="none", # trace lines inside the heatmap
                           cexRow=1.1,
                           cexCol=1.1,
                           col=map.colour,
                           #breaks=palette.breaks, # pairs.breaks or "none"

                           margins = c(10,24), # extra Y and X margins respectively # 9 and 18
                           #lhei = c(1,4), # lhei = c(1,7), # relative height of the legend and plot respectively # default is c(1.5,4,1)
                           #lwid = c(1,4), # lwid = c(1,4), # relativewidth of the legend and plot respectively # default is c(1.5,4)

                           sepwidth = c(0.1, 0.1),
                           sepcolor = "white",
                           #rowsep = c(6, 8, 10, 14, 17, 20, 21, 23, 24),
                           #colsep = c(1, 9, 10, 14, 20, 22, 24, 28, 31),

                           density.info="histogram") # "histogram" "density" "none"

      ## CLICK EXPORT / CHOOSE DIRECTORY TO SAVE / SAVE AS PDF
      ## Choose A4 and LANDSCAPE -- re-export with modified dimensions if required

    }








make.heatmap.fold <- function(x,
                              sample.col,
                              col.rmv
                              plot.title,
                              do.transpose,
                              max.range,
                              mix.range,
                              do.dendrogram,
                              dendro.set,
                              row.clustering
                              n.row.groups,
                              col.clustering,
                              n.col.groups)





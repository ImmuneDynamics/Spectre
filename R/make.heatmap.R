#' make.heatmap - Create a heatmap
#'
#' This function allows you to create a coloured heatmap for clusters vs markers (coloured by MFI) or clusters vs samples (coloured by MFI or cell counts). Can be set to plot 'fold-change' values that are provided in log2.
#' make.heatmap is a wrapper around the heatmap.2 function from gplots.
#'
#' @usage make.heatmap(x, ...)
#'
#' @param x data.frame. Clusters/populations vs markers (MFI) or Clusters/populations vs samples or (MFI or cell numbers). No default.
#' @param sample.col Character. Specify the name of the column that indicates samples. No default.
#' @param annot.cols Character. Columns which contain values that you do NOT want to plot in the heatmap, e.g. sample name, group name, Live/Dead etc. No default.
#' @param plot.title Character.
#' @param transpose Logical. Do you want to transpose the heatmap. Defaults to FALSE.
#' @param is.fold Logical. TRUE for fold-change heatmap, FALSE for normal values heatmap. Defaults to FALSE.
#' @param fold.max.range Numeric. For fold-change heatmaps, what is the maxium colour value that should be plotted. Defaults to 3.
#' @param fold.min.range Numeric. For fold-change heatmaps, what is the minimum colour value that should be plotted. Defaults to -3.
#' @param normalise Logical. Only applies to standard heatmaps (i.e. when is.fold = FALSE). TRUE to normalise each column between 0 and 1, FALSE to plot the raw values. Defaults to TRUE.
#' @param dendrograms Character. Do you want to include dendrograms on columns/rows. Can be "both", "row", "column", or "none. Defaults to "both.
#' @param n.row.groups Numeric.
#' @param n.col.groups Numeric.
#' @param row.sep Numeric.
#' @param col.sep Numeric.
#' @param y.margin Numeric.
#' @param x.margin Numeric.
#' @param standard.colours Character.
#' @param fold.colours Character.
#' @param colour.scheme Character.
#'
#' @return prints a heatmap with dendrograms.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @examples
#' make.heatmap()
#'
#' @export

# to add
  # extra margins for column and row names
  # calculation -- recommended W and H based on how man rows and columns there are -- to keep each heatmap cell a square
  #

make.heatmap <- function(x,
                         sample.col,
                         annot.cols,
                         plot.title,

                         transpose        = FALSE,

                         is.fold          = FALSE,
                         fold.max.range    = 3,
                         fold.min.range    = -3,

                         normalise        = TRUE,

                         dendrograms      = "both",
                         n.row.groups     = 2,   # choose number of row groupings (by default, rows are samples, unless do.transpose = 1). Must be >1.
                         n.col.groups     = 2,   # choose number of column groupings (by default, columns are populations, unless do.transpose = 1). Must be >1.

                         row.sep          = c(),
                         col.sep          = c(),

                         y.margin         = 1,
                         x.margin         = 1,

                         standard.colours = "colour.palette", # default and only current option
                         fold.colours     = "fold.palette",   # default and only current option
                         colour.scheme    = "group.palette")  # default and only current option
                         #plot.width       = 11.69,    # Default = 11.69 inches, target for A4 landscape. # Default = 8.26 inches, target for A4 landscape.
                         #plot.height      = 8.26)
{

  ### TESTING

      # library(gplots)
      #
      # setwd("/Users/Tom/Google Drive (t.ashhurst@centenary.org.au)/_Sydney Cytometry/2019_Synced/GitHub/Public github/Spectre/Workflow scripts/Output_CAPX/Output-SumTables/Output_CellNums/")
      #
      # x <- read.csv("SumTable_CellsPerTissue.csv")
      # sample.col <- "Sample"
      # annot.cols <- c(1:3)
      # plot.title <- "Test"
      #
      # transpose        = FALSE
      # is.fold          = FALSE
      # normalise        = TRUE
      # dendrograms      = "row"
      #
      # n.row.groups     = 2
      # n.col.groups     = 2
      #
      # fold.max.range    = 3
      # fold.min.range    = -3
      #
      # y.margin = 8
      # x.margin = 8
      #
      # row.sep          = c(6)
      # col.sep          = c(6)
      # standard.colours = "colour.palette"
      # fold.colours     = "fold.palette"
      # colour.scheme    = "group.palette"
      #
      # plot.width       = 11.69
      # plot.height      = 8.26

  ### Setup color scheme

      ## Standard colour options
      colour.palette <- (colorRampPalette(brewer.pal(9, "YlGnBu"))(31)) # 256

      ## Fold-change colour options
      fold.palette <- colorRampPalette(rev(c("#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026","black","#023858","#045a8d","#0570b0","#3690c0","#74a9cf","#a6bddb","#d0d1e6","#ece7f2")))

      ## Grouping colour options
      group.palette <- (colorRampPalette(brewer.pal(11, "Spectral")))
      white.palette <- colorRampPalette(c("white"))

  ### 2.1 - Embed cluster or population name as row name

      heatmap.data <- x
      rownames(heatmap.data) <- t(x[sample.col])
      heatmap.data

      heatmap.data <- heatmap.data[-c(annot.cols)] # remove columns
      heatmap.data

  ### 2.2 - Transpose (ONLY IF REQUIRED) -- the longest set (clusters or parameters) on x-axis -- by default MARKERS are columns, CLUSTERS are rows -- transpose to flip these defaults
      if(transpose == TRUE){
        heatmap.data.t <- as.data.frame(t(heatmap.data))
        heatmap.data <- heatmap.data.t
      }

  ### 2.3 - NORMALISE BY COLUMN (i.e. each column/parameter has a max of 1 and a minimum of 0) # This is optional, but allows for better comparison between markers
      if(normalise == TRUE){
        if(is.fold == FALSE){
          row.nam <- row.names(heatmap.data)

          norm.fun <- function(x) {(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) -min(x, na.rm=TRUE))}
          heatmap.data.norm <- as.data.frame(lapply(heatmap.data, norm.fun)) # by default, removes the names of each row
          max(heatmap.data.norm)
          max(heatmap.data.norm[,5])
          heatmap.data.norm <- as.matrix(heatmap.data.norm)
          heatmap.data <- heatmap.data.norm
          rownames(heatmap.data) <- row.nam # add row names back
        }
      }

      # convert to matrix
      heatmap.data <- as.matrix(heatmap.data)

  ### 2.4 - Set up clustering

      if(dendrograms == "none"){
        row.clustering    <- FALSE
        col.clustering    <- FALSE
      }

      if(dendrograms != "none"){

        # set the custom distance and clustering functions, per your example
        hclustfunc <- function(x) hclust(x, method="complete")
        distfunc <- function(x) dist(x, method="euclidean")

        # perform clustering on rows and columns
            if(dendrograms == "both"){
              row.clustering    <- TRUE
              col.clustering    <- TRUE

              cl.row <- hclustfunc(distfunc(heatmap.data))
              cl.col <- hclustfunc(distfunc(t(heatmap.data)))

              nrow(heatmap.data)# work out no cols and rows
              ncol(heatmap.data)

              gr.row <- cutree(cl.row, n.row.groups) # extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
              gr.col <- cutree(cl.col, n.col.groups)

              col1 <- group.palette(n.row.groups)
              col2 <- group.palette(n.col.groups)
            }

            if(dendrograms == "column"){
              row.clustering    <- FALSE
              col.clustering    <- TRUE

              cl.row <- hclustfunc(distfunc(heatmap.data))
              cl.col <- hclustfunc(distfunc(t(heatmap.data)))

              nrow(heatmap.data)# work out no cols and rows
              ncol(heatmap.data)

              gr.row <- cutree(cl.row, n.row.groups) # extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
              gr.col <- cutree(cl.col, n.col.groups)

              col1 <- white.palette(n.row.groups)
              col2 <- group.palette(n.col.groups)
              }

            if(dendrograms == "row"){
              row.clustering    <- TRUE
              col.clustering    <- FALSE

              cl.row <- hclustfunc(distfunc(heatmap.data))
              cl.col <- hclustfunc(distfunc(t(heatmap.data)))

              nrow(heatmap.data)# work out no cols and rows
              ncol(heatmap.data)

              gr.row <- cutree(cl.row, n.row.groups) # extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
              gr.col <- cutree(cl.col, n.col.groups)

              col1 <- group.palette(n.row.groups)
              col2 <- white.palette(n.col.groups)
            }
      }

  ### 2.6 - Set up fold-change or normal preferences

      if(is.fold == TRUE){
        map.colour <- fold.palette(31)
        sym.key <- FALSE # TRUE need if NOT creating own breaks
        sym.breaks <- TRUE
        my.breaks <- seq(fold.min.range, fold.max.range, length.out = 32)
      }

      if(is.fold == FALSE){
        map.colour <- colour.palette
        sym.key <- FALSE
        sym.breaks <- FALSE
        my.breaks <- seq(min(heatmap.data), max(heatmap.data), length.out = 32)
      }

      scale.set <- "none" # Can be: "none", "both", "row", "column"

      #if(is.fold == TRUE){a <- "(fold-change)"}
      #if(is.fold == FALSE){a <- ""}

      #title.text <- paste0(plot.title, a, ".pdf")
      title.text <- plot.title

  ### 2.7 - Plot heatmap and dendrograms

      ## dev.off()
      par(cex.main=1.3)
      par(xpd = TRUE)

      # pdf(title.text, width = plot.width, height = plot.height)
      # HeatMap <-
      print(heatmap.2(as.matrix(heatmap.data),  #### http://stanstrup.github.io/heatmaps.html # https://stackoverflow.com/questions/22278508/how-to-add-colsidecolors-on-heatmap-2-after-performing-bi-clustering-row-and-co
                           main=title.text,
                           notecol= "black", # font colour of all cells to black
                           key=TRUE,
                           #keysize=0.75,#key.title = NA,
                           #key.par = list(cex=1),

                           hclust=hclustfunc,
                           distfun=distfunc,

                           dendrogram = dendrograms, # "both", "row", "column" etc
                           Rowv = row.clustering, # turns on row clustering # = as.dendrogram(cluster)
                           Colv = col.clustering, # turns on column clustering # = reorder(as.dendrogram(cluster.row), 10:1) to flip around, or maybe dictage

                           RowSideColors = col1[gr.row],
                           ColSideColors = col2[gr.col], # as.character(clusters)

                           revC = FALSE, # default FALSE
                           symm = FALSE, # default FALSE
                           symkey = sym.key, # default FALSE, TRUE FOR SYMMETRY IN FOLD CHANGE SITUATIONS
                           symbreaks = sym.breaks, # default FALSE

                           breaks = my.breaks,

                           scale = scale.set, # "none" or "row" or "column". default is to scale by row, I think? # scale argument only for heat data, not dendrogram
                           trace = "none", # trace lines inside the heatmap
                           cexRow = 1.1,
                           cexCol = 1.1,
                           col = map.colour,
                           #breaks=palette.breaks, # pairs.breaks or "none"

                           #margins = c(9,18), # extra Y and X margins respectively # 9 and 18
                           #lhei = c(1,4), # lhei = c(1,7), # relative height of the legend and plot respectively # default is c(1.5,4,1)
                           #lwid = c(1,4), # lwid = c(1,4), # relativewidth of the legend and plot respectively # default is c(1.5,4)

                           margins = c(y.margin, x.margin), # extra Y and X margins respectively # 9 and 18

                           sepwidth = c(0.1, 0.1),
                           sepcolor = "white",
                           ##### FIX
                           rowsep = row.sep,
                           colsep = col.sep,
                           #####

                           density.info="histogram")) # "histogram" "density" "none"
      #dev.off()

}

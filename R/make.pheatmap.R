#' make.pheatmap - Create a 'pretty' heatmap (pheatmap)
#'
#' This function allows you to create a coloured heatmap for clusters vs markers (coloured by MFI) or clusters vs samples (coloured by MFI or cell counts).
#' Can be set to plot 'fold-change' values (with is.fold = TRUE) that are provided in log2.
#' make.pheatmap is a wrapper around the pheatmap function from pheatmap.
#' Also uses the packages 'RColorBrewer' and 'scales' for colour customisation.
#'
#' @param dat NO DEFAULT. data.frame. Clusters/populations vs markers (MFI) or Clusters/populations vs samples or (MFI or cell numbers).
#' @param sample.col NO DEFAULT. Character. Specify the name of the column that indicates samples. Will be rows on heatmap.
#' @param plot.cols NO DEFAULT. Character vector. Name of columns you wish to plot on the heatmap. Will be columns on heatmap.
#' @param annot.cols DEFAULT = NULL. Character. Columns which contain values that you do NOT want to plot in the heatmap, e.g. sample name, group name, Live/Dead etc.
#' @param file.name DEFAULT = paste0("Pheatmap by ", sample.col, ".png"). Character. What do you want to call the file, including the extension.
#' @param plot.title DEFAULT = paste0(sample.col, " heatmap"). Character.
#' @param transpose DEFAULT = FALSE. Logical. Do you want to transpose the heatmap.
#' @param is.fold DEFAULT = FALSE. Logical. TRUE for fold-change heatmap, FALSE for normal values heatmap.
#' @param fold.range DEFAULT = NULL. Numeric vector. For fold-change heatmaps, what is the maxium and minimum values that should be plotted. Example: for a max of 3, and a minimum of -3 would b c(3,-3). Defaults to NULL (which will use the max and min within the dataset).
#' @param normalise DEFAULT = TRUE. Logical. Only applies to standard heatmaps (i.e. when is.fold = FALSE). TRUE to normalise each column between 0 and 1, FALSE to plot the raw values.
#' @param dendrograms DEFAULT = "both". Character. Do you want to include dendrograms on columns/rows. Can be "both", "row", "column", or "none.
#' @param cutree_rows DEFAULT = 1 (i.e. no divisions). Divides rows based on the level of dendrogram branches, if dendrograms = 'both' or 'row'.
#' @param cutree_cols DEFAULT = 1 (i.e. no divisions). Divides columns based on the level of dendrogram branches, if dendrograms = 'both' or 'column'
#' @param row.sep DEFAULT = c(). Numeric. Only used if not clustering rows.
#' @param col.sep DEFAULT = c(). Numeric. Only used if not clustering columns
#' @param cell.size DEFAULT = 15. Numeric.
#' @param standard.colours DEFAULT = "BuPu". Character. Can also be "RdYlBu", "YlGnBu", "viridis", "magma", "inferno", "spectral", "Blues", "Reds", "Greys", or "rev(RdBu)".
#' @param fold.colours DEFAULT = "Spectre". Character, a custom blue:black:red fold-change colour scheme. Can also be "RdYlBu", "YlGnBu", "viridis", "magma", "inferno", "spectral", "Blues", "Reds", "Greys", or "rev(RdBu)".
#' @param path DEFAULT = NULL. The location to save plots. By default, will save to current working directory. Can be overidden by entering a sub-directory name (e.g. if you have a subdirectory under your working directory called 'plots', you can enter this here).
#'
#' @usage make.pheatmap(dat, sample.col, plot.cols, annot.cols, file.name, plot.title, transpose, is.fold, fold.range, normalise, dendrograms, row.sep, col.sep, cell.size, standard.colours, path)
#'
#' @examples
#' ## MFI cluster vs marker heatmap
#' Spectre::make.pheatmap(dat = Spectre::demo.exp,
#'                        file.name = "Expression pheatmap.png",
#'                        plot.title = "Expression",
#'                        sample.col = "Population",
#'                        plot.cols = names(Spectre::demo.exp)[c(2:10)])
#'
#' ## Z-scrore of fold-change type heatmap
#' z.dat <- do.zscore(dat = Spectre::demo.sum,
#'                    use.cols = names(Spectre::demo.sum)[c(4:15)],
#'                    replace = TRUE)
#'
#' Spectre::make.pheatmap(dat = z.dat,
#'                        file.name = "z-score.png",
#'                        plot.title = "z-score",
#'                        sample.col = "Sample",
#'                        plot.cols = names(z.dat)[c(4:15)],
#'                        annot.cols = names(z.dat)[c(2:3)],
#'                        is.fold = TRUE,
#'                        fold.range = c(3,-3)
#'                        )
#'
#' @author
#' Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#'
#' @references
#' \url{https://sydneycytometry.org.au/spectre}.
#' Helpful examples at \url{https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/}
#'
#' @import data.table
#'
#' @export

make.pheatmap <- function(dat,
                          
                          sample.col,
                          plot.cols,
                          
                          annot.cols = NULL,
                          feature.annots = NULL,
                          annotation_colors = NULL,
                          
                          file.name = paste0("Pheatmap by ", sample.col, ".png"),
                          plot.title = paste0(sample.col, " heatmap"),
                          
                          transpose = FALSE,
                          normalise = TRUE,
                          is.fold = FALSE,
                          fold.range = NULL, # c(3,-3)
                          dendrograms = "both",
                          
                          dendrograms.sort = FALSE, # if TRUE, will SORT dendrograms
                          
                          cutree_rows = 1,
                          cutree_cols = 1,
                          row.sep = c(),
                          col.sep = c(),
                          cell.size = 15,
                          standard.colours = "BuPu",
                          fold.colours = 'Spectre',
                          path = NULL)

#n.row.groups = 2,
#n.col.groups = 2,
#plot.width       = 11.69,    # Default = 11.69 inches, target for A4 landscape. # Default = 8.26 inches, target for A4 landscape.
#plot.height      = 8.26)
{
  
  ### Check that necessary packages are installed
  if(!is.element('pheatmap', installed.packages()[,1])) stop('pheatmap is required but not installed')
  if(!is.element('RColorBrewer', installed.packages()[,1])) stop('RColorBrewer is required but not installed')
  if(!is.element('scales', installed.packages()[,1])) stop('scales is required but not installed')
  
  ### Require packages
  require(pheatmap)
  require(RColorBrewer)
  require(scales)
  
  ### Setup color scheme
  
  ## Standard colour options
  
  if(standard.colours == "BuPu"){
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, "BuPu"))(31)) # 256
  }
  
  if(standard.colours == "RdYlBu"){
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu"))(31)) # 256
    colour.palette <- rev(colour.palette)
  }
  
  if(standard.colours == "rev(RdBu)"){
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(31)) # 256
    colour.palette <- rev(colour.palette)
  }
  
  if(standard.colours == "Blues"){
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(31)) # 256
  }
  
  if(standard.colours == "Reds"){
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))(31)) # 256
  }
  
  if(standard.colours == "Greys"){
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, "Greys"))(31)) # 256
  }
  
  if(standard.colours == "YlGnBu"){
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(31)) # 256
  }
  
  if(standard.colours == "viridis"){
    colour.palette <- colorRampPalette(c(scales::viridis_pal(option = "viridis")(50)))
    colour.palette <- colour.palette(31)
  }
  
  if(standard.colours == "spectral"){
    spectral.list <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(50)
    spectral.list <- rev(spectral.list)
    colour.palette <- colorRampPalette(c(spectral.list))
    colour.palette <- colour.palette(31)
  }
  
  if(standard.colours == "magma"){
    colour.palette <- colorRampPalette(c(scales::viridis_pal(option = "magma")(50)))
    colour.palette <- colour.palette(31)
  }
  
  if(standard.colours == "inferno"){
    colour.palette <- colorRampPalette(c(scales::viridis_pal(option = "inferno")(50)))
    colour.palette <- colour.palette(31)
  }
  
  ## Fold-change colour options
  
  if(fold.colours == "Spectre"){
    fold.palette <- colorRampPalette(rev(c("#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026","black","#023858","#045a8d","#0570b0","#3690c0","#74a9cf","#a6bddb","#d0d1e6","#ece7f2")))
    fold.palette <- fold.palette(31)
  }
  
  if(fold.colours == "BuPu"){
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, "BuPu"))(31)) # 256
  }
  
  if(fold.colours == "RdYlBu"){
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu"))(31)) # 256
    fold.palette <- rev(fold.palette)
  }
  
  if(fold.colours == "rev(RdBu)"){
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(31)) # 256
    fold.palette <- rev(fold.palette)
  }
  
  if(fold.colours == "Blues"){
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(31)) # 256
  }
  
  if(fold.colours == "Reds"){
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))(31)) # 256
  }
  
  if(fold.colours == "Greys"){
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, "Greys"))(31)) # 256
  }
  
  if(fold.colours == "YlGnBu"){
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(31)) # 256
  }
  
  if(fold.colours == "viridis"){
    fold.palette <- colorRampPalette(c(scales::viridis_pal(option = "viridis")(50)))
    fold.palette <- fold.palette(31)
  }
  
  if(fold.colours == "spectral"){
    spectral.list <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(50)
    spectral.list <- rev(spectral.list)
    fold.palette <- colorRampPalette(c(spectral.list))
    fold.palette <- fold.palette(31)
  }
  
  if(fold.colours == "magma"){
    fold.palette <- colorRampPalette(c(scales::viridis_pal(option = "magma")(50)))
    fold.palette <- fold.palette(31)
  }
  
  if(fold.colours == "inferno"){
    fold.palette <- colorRampPalette(c(scales::viridis_pal(option = "inferno")(50)))
    fold.palette <- fold.palette(31)
  }
  
  #fold.palette <- colorRampPalette(brewer.pal(11, "RdBu"))
  
  ## Grouping colour options
  # group.palette <- (colorRampPalette(brewer.pal(11, "Spectral")))
  # white.palette <- colorRampPalette(c("white"))
  
  ### Embed cluster or population name as row name
  
  dat <- as.data.frame(dat)
  
  heatmap.data <- dat
  rownames(heatmap.data) <- t(dat[sample.col])
  heatmap.data
  
  if(is.null(annot.cols) == FALSE){
    annot <- heatmap.data[annot.cols]
    heatmap.data <- heatmap.data[plot.cols] # remove columns
    heatmap.data
  }
  
  if(is.null(annot.cols) == TRUE){
    annot <- NULL
    heatmap.data <- heatmap.data[plot.cols] # remove columns
    heatmap.data
  }
  
  ### Transpose (ONLY IF REQUIRED) -- the longest set (clusters or parameters) on x-axis -- by default MARKERS are columns, CLUSTERS are rows -- transpose to flip these defaults
  if(transpose == TRUE){
    heatmap.data.t <- as.data.frame(t(heatmap.data))
    heatmap.data <- heatmap.data.t
  }
  
  ### NORMALISE BY COLUMN (i.e. each column/parameter has a max of 1 and a minimum of 0) # This is optional, but allows for better comparison between markers
  if(normalise == TRUE){
    if(is.fold == FALSE){
      row.nam <- row.names(heatmap.data)
      col.nam <- names(heatmap.data)
      
      norm.fun <- function(x) {(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) -min(x, na.rm=TRUE))}
      heatmap.data.norm <- as.data.frame(lapply(heatmap.data, norm.fun)) # by default, removes the names of each row
      names(heatmap.data.norm) <- col.nam
      
      max(heatmap.data.norm)
      # max(heatmap.data.norm[,5])
      heatmap.data.norm <- as.matrix(heatmap.data.norm)
      heatmap.data <- heatmap.data.norm
      rownames(heatmap.data) <- row.nam # add row names back
    }
  }
  
  # convert to matrix
  heatmap.data <- as.matrix(heatmap.data)
  
  ### Set up clustering
  
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
      
      if(isTRUE(dendrograms.sort)){
        
        row.clustering <- hclustfunc(distfunc(heatmap.data))
        col.clustering <- hclustfunc(distfunc(t(heatmap.data)))
        
        require(dendsort)
        sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
        
        row.clustering <- sort_hclust(row.clustering)
        col.clustering <- sort_hclust(col.clustering)
      }
      
      # cl.row <- hclustfunc(distfunc(heatmap.data))
      # cl.col <- hclustfunc(distfunc(t(heatmap.data)))
      
      # nrow(heatmap.data)# work out no cols and rows
      # ncol(heatmap.data)
      
      # gr.row <- cutree(cl.row, n.row.groups) # extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
      # gr.col <- cutree(cl.col, n.col.groups)
      
      # col1 <- group.palette(n.row.groups)
      # col2 <- group.palette(n.col.groups)
    }
    
    if(dendrograms == "column"){
      row.clustering    <- FALSE
      col.clustering    <- TRUE
      
      if(isTRUE(dendrograms.sort)){
        
        col.clustering <- hclustfunc(distfunc(t(heatmap.data)))
        
        require(dendsort)
        sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
        
        col.clustering <- sort_hclust(col.clustering)
      }
      
      # cl.row <- hclustfunc(distfunc(heatmap.data))
      # cl.col <- hclustfunc(distfunc(t(heatmap.data)))
      
      # nrow(heatmap.data)# work out no cols and rows
      # ncol(heatmap.data)
      
      # gr.row <- cutree(cl.row, n.row.groups) # extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
      # gr.col <- cutree(cl.col, n.col.groups)
      
      # col1 <- white.palette(n.row.groups)
      # col2 <- group.palette(n.col.groups)
    }
    
    if(dendrograms == "row"){
      row.clustering    <- TRUE
      col.clustering    <- FALSE
      
      if(isTRUE(dendrograms.sort)){
        
        row.clustering <- hclustfunc(distfunc(t(heatmap.data)))
        
        require(dendsort)
        sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
        
        row.clustering <- sort_hclust(row.clustering)
      }
      
      # cl.row <- hclustfunc(distfunc(heatmap.data))
      # cl.col <- hclustfunc(distfunc(t(heatmap.data)))
      
      # nrow(heatmap.data)# work out no cols and rows
      # ncol(heatmap.data)
      
      # gr.row <- cutree(cl.row, n.row.groups) # extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
      # gr.col <- cutree(cl.col, n.col.groups)
      
      # col1 <- group.palette(n.row.groups)
      # col2 <- white.palette(n.col.groups)
    }
  }
  
  ### Set up fold-change or normal preferences
  
  if(is.fold == TRUE){
    map.colour <- fold.palette
    sym.key <- FALSE # TRUE need if NOT creating own breaks
    sym.breaks <- TRUE
    
    if(is.null(fold.range)){
      fld.max <- max(heatmap.data, na.rm = TRUE)
      fld.min <- min(heatmap.data, na.rm = TRUE)
      
      if(fld.max == -fld.min){
        fold.max.range <- fld.max
        fold.min.range <- fld.min
      }
      
      if(fld.max > -fld.min){
        fold.max.range <- fld.max
        fold.min.range <- -fld.max
      }
      
      if(fld.max < -fld.min){
        fold.max.range <- -fld.min
        fold.min.range <- fld.min
      }
    }
    
    if(!is.null(fold.range)){
      fold.max.range <- fold.range[1]
      fold.min.range <- fold.range[2]
    }
    
    my.breaks <- seq(fold.min.range, fold.max.range, length.out = 32)
    
  }
  
  if(is.fold == FALSE){
    map.colour <- colour.palette
    sym.key <- FALSE
    sym.breaks <- FALSE
    heatmap.data
    
    my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
    my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
    
    my.breaks <- seq(my.min(heatmap.data), my.max(heatmap.data), length.out = 32)
  }
  
  scale.set <- "none" # Can be: "none", "both", "row", "column"
  
  #if(is.fold == TRUE){a <- "(fold-change)"}
  #if(is.fold == FALSE){a <- ""}
  
  #title.text <- paste0(plot.title, a, ".pdf")
  title.text <- plot.title
  
  ### Plot heatmap and dendrograms
  
  # Specify directory heatmap will be saved
  
  if(is.null(path)){
    flnm <- file.name
  }
  
  if(!is.null(path)){
    flnm <- paste0(path, '/', file.name)
  }
  
  pheatmap::pheatmap(mat = as.matrix(heatmap.data),
                     main = title.text,
                     
                     cellwidth = cell.size,
                     cellheight = cell.size,
                     
                     cluster_rows = row.clustering,
                     cluster_cols = col.clustering,
                     #scale = "column",
                     
                     breaks = my.breaks,
                     
                     cutree_rows = cutree_rows,
                     cutree_cols = cutree_cols,
                     
                     gaps_row = row.sep,
                     gaps_col = col.sep,
                     
                     annotation_row = annot,
                     annotation_col = feature.annots,
                     annotation_colors = annotation_colors,
                     
                     color = map.colour,
                     filename = flnm
  )
  
  message(paste0("A pheatmap has been saved to your working directory", paste0(path, file.name)))
  
}




# if(!require('pheatmap')) {install.packages('pheatmap')}
# library(pheatmap)
#
# size <- 10
#
# pheatmap(mat = B, #A[-c(1:3)]
#          cellwidth = size,
#          cellheight = size,
#
#          cluster_rows = TRUE,
#          cluster_cols = FALSE,
#
#          #scale = "column",
#          main = "Test Heatmap",
#          #color = rev(c("#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026","black","#023858","#045a8d","#0570b0","#3690c0","#74a9cf","#a6bddb","#d0d1e6","#ece7f2")),
#          filename = "testPheatmap.png"
#          )
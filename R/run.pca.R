#' Run the PCA algorithm (using stats::prcomp)
#' 
#' Method to run a PCA dimensionality reduction algorithm.
#' A principal component analysis (PCA) is capable of reducing the number of dimensions (i.e. parameters) with minimal effect on the variation of the given dataset.
#' This function will run a PCA calculation (extremely fast) and generate plots (takes time).
#' For individuals (such as samples or patients), a PCA can group them based on their similarities.
#' A PCA is also capable of ranking variables/parameters (such as markers or cell counts) based on their contribution to the variability across a dataset in an extremely fast manner.
#' In cytometry, this can be useful to identify marker(s) that can be used to differentiate between subset(s) of cells.
#' Uses the base R package "stats" for PCA, "factoextra" for PCA and scree plots, "data.table" for saving .csv files, "ggplot2" for saving plots, "gtools" for rearranging data order, 'RColorBrewer' and 'viridis' for colour schemes.
#' More information on PCA plots can be found here http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/.
#'
#' @param dat NO DEFAULT. data.frame.
#' @param use.cols NO DEFAULT. Vector of numbers, reflecting the columns to use for dimensionality reduction (may not want parameters such as "Time" or "Sample").
#' @param scale DEFAULT = TRUE. A logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place.
#' @param scree.plot DEFAULT = TRUE. Option to create scree plots. Note this will require the input of an elbow point during run. Will save generated scree plot.
#' @param variable.contribution DEFAULT = TRUE. Option to create plot showing the contribution of each variable. Horizontal red line represents the average variable contribution if all variables contributed equally. Requires scree.plot = TRUE.
#' @param plot.individuals DEFAULT = TRUE. Option to create PCA plots on individuals (samples/patients).
#' @param plot.ind.label DEFAULT = "point". Option to add text to PCA plots on individuals as an extra identifier. Use c("point", "text") to include both text and point.
#' @param pointsize.ind DEFAULT = 1.5. Numeric. Size of dots of individuals on PCA plot.
#' @param row.names DEFAULT = NULL. Column (as character) that defines individuals. Will be used to place name on plot.individuals.
#' @param plot.ind.group DEFAULT = FALSE. Option to group inidividuals with ellipses (which by default show the 95 % confidence interval). Must specify column that groups individuals with group.ind.
#' @param group.ind DEFAULT = NULL. Column (as character) that defines groups of individuals. Works with plot.ind.group which must be set to TRUE.
#' @param colour.group DEFAULT = "viridis". Colour scheme for each group. Options include "jet", "spectral", "viridis", "inferno", "magma".
#' @param pointsize.group DEFAULT = 1.5. Numeric. Size of shapes of group individuals on PCA plot.
#' @param ellipse.type DEFAULT = "confidence". Set type of ellipse. Options include "confidence", "convex", "concentration", "t", "norm", "euclid". See factoextra::fviz for more information.
#' @param ellipse.level DEFAULT = 0.95. Size of ellipses. By default 95 % (0.95).
#' @param plot.variables DEFAULT = TRUE. Option to create PCA plots on variables (markers/cell counts).
#' @param colour.var DEFAULT = "solid". Colour scheme for PCA plot with variables. Options include "solid", "jet", "spectral", "viridis", "inferno", "magma", "BuPu". Note some colours are pale and may not appear clearly on plot.
#' @param plot.combined DEFAULT = TRUE. Option to create a combined PCA plot with both individuals and variables.
#' @param repel DEFAULT = FALSE. Option to avoid overlapping text in PCA plots. Can greatly increase plot time if there is a large number of samples.
#' @param var.numb DEFAULT = 20. Top number of variables to be plotted. Note the greater the number, the longer plots will take.
#' @param path DEFAULT = getwd(). The location to save plots. By default, will save to current working directory. Can be overidden.
#' 
#' @usage run.pca(dat, use.cols, scale = TRUE, scree.plot = TRUE, variable.contribution = TRUE, plot.individuals = TRUE, plot.ind.label = "point", pointsize.ind = 1.5, row.names = NULL, plot.ind.group = FALSE, group.ind = NULL, colour.group = "viridis", pointsize.group = 1.5, ellipse.type = "confidence", ellipse.level = 0.95, plot.variables = TRUE, colour.var = "solid", plot.combined = TRUE, repel = FALSE, var.numb = 20, path = getwd())
#'
#' @examples
#' # Set directory to save files. By default it will save files at get()
#' setwd("/Users/felixmarsh-wakefield/Desktop")
#' 
#' # Run PCA on demonstration dataset
#' Spectre::run.pca(dat = Spectre::demo.clustered,
#'                 use.cols = c(11:19),
#'                 repel = TRUE
#'                 )
#' 
#' # Compare between groups
#' Spectre::run.pca(dat = Spectre::demo.clustered,
#'                  use.cols = c(11:19),
#'                  plot.ind.label = c("point", "text"), #individual cells will be labelled as numbers
#'                  plot.ind.group = TRUE,
#'                  group.ind = "Group"
#'                  )
#'         
#' # When prompted, type in "5" and click enter to continue function (this selects the elbow point based off the scree plot)
#' 
#' ## Possible issues ##
#' # Remove any NA present
#' na.omit(dat)
#' 
#' # Remove columns that have zero variance (e.g. if MFI is the same for all samples for a marker)
#' dat <- data.table::as.data.table(dat)
#' dat <- dat[ , lapply(.SD, function(v) if(data.table::uniqueN(v, na.rm = TRUE) > 1) v)] #for data table format
#' 
#' # Ellipses are only generated in 'plot.ind.group' when there are at least 2 samples per group ('group.ind')
#' 
#' @author Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#' @export

run.pca <- function(dat,
                   use.cols,
                   scale = TRUE,
                   scree.plot = TRUE,
                   variable.contribution = TRUE,
                   plot.individuals = TRUE,
                   plot.ind.label = "point",
                   pointsize.ind = 1.5,
                   row.names = NULL,
                   plot.ind.group = FALSE,
                   group.ind = NULL,
                   colour.group = "viridis",
                   pointsize.group = 1.5,
                   ellipse.type = "confidence",
                   ellipse.level = 0.95,
                   plot.variables = TRUE,
                   colour.var  = "solid",
                   plot.combined = TRUE,
                   repel = FALSE,
                   var.numb = 20,
                   path = getwd()
                   ) {
  
  ## Check that necessary packages are installed
  if(!is.element('stats', installed.packages()[,1])) stop('stats is required but not installed')
  if(!is.element('factoextra', installed.packages()[,1])) stop('factoextra is required but not installed')
  if(!is.element('ggplot2', installed.packages()[,1])) stop('ggplot2 is required but not installed')
  if(!is.element('gtools', installed.packages()[,1])) stop('gtools is required but not installed')
  if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')
  if(!is.element('RColorBrewer', installed.packages()[,1])) stop('RColorBrewer is required but not installed')
  if(!is.element('viridis', installed.packages()[,1])) stop('viridis is required but not installed')
  
  ## Require packages
  require(stats)
  require(factoextra)
  require(ggplot2)
  require(gtools)
  require(data.table)
  require(RColorBrewer)
  require(viridis)
  
  # Make sure dat is a dataframe
  dat <- as.data.frame(dat)
  
  # Change rownames of dat
  if (is.null(row.names) == FALSE) {
    rownames(dat) <- dat[, row.names]
  }
  
  ## Run PCA
  pca_out <- stats::prcomp(dat[, use.cols],
                           scale = scale
                           )
  
  # Create scree plot
  if (scree.plot == TRUE) {
    scree_plot <- factoextra::fviz_eig(pca_out, addlabels = TRUE) #creates scree plot; addlabels adds % to plot
    print(scree_plot)
    
    elbow.point <- readline("Type in the elbow point based on the scree plot. Must be positive integer. ")
    elbow.point <- as.numeric(elbow.point)
    
    # Saves scree plot
    ggplot2::ggsave(scree_plot,
           filename = "Scree plot.pdf",
           path = path,
           useDingbats = FALSE)
    
    if (variable.contribution == TRUE) {
      # Creates plot with contributions of each variable
      pca.contrib <- factoextra::fviz_contrib(pca_out, #red dashed line represents the expected average row contributions if the contributions were uniform: 1/nrow(markers)
                                  choice = "var", #variable ("var") or individual ("ind")
                                  axes = 1:elbow.point, #axes refer to PCA number
                                  top = var.numb #top number of variables to be plotted
                                  )
      
      ggplot2::ggsave(pca.contrib,
             filename = "PCA plot-contribution.pdf",
             width = 16,
             height = 5,
             path = path,
             useDingbats = FALSE)
      
      # Extract contribution of variables
      pca_out_var <- factoextra::get_pca_var(pca_out)
      pca_var_contrib <- as.data.frame(pca_out_var$contrib) #coordinates of each point for PCA plot

      data.table::fwrite(pca_var_contrib,
                         "PCA-contributions.csv",
                         row.names = TRUE)
      
      # Extract raw eigenvalues
      pca_eig <- factoextra::get_eig(pca_out)
      
      data.table::fwrite(pca_eig,
                         "PCA-eigenvalues.csv",
                         row.names = TRUE)
      
      # Extract corrected contributions
      pca_eig_contrib <- facto_summarize(pca_out,
                                         element = "var",
                                         result = "contrib",
                                         axes = 1:elbow.point
                                         )
      
      data.table::fwrite(pca_eig_contrib,
                         paste0("PCA-eig-contrib_", elbow.point, "-dim.csv"),
                         row.names = TRUE)
      
    }
    
  }
  
  # Create PCA plot with individuals
  if (plot.individuals == TRUE) {
    pca_plot_ind <- factoextra::fviz_pca_ind(pca_out,
                                             repel = repel,
                                             geom = plot.ind.label,
                                             pointsize = pointsize.ind
                                             )
    
    ggplot2::ggsave(pca_plot_ind,
                    filename = "PCA plot-individuals.pdf",
                    path = path,
                    useDingbats = FALSE)
    
    # Extract coordinates of individuals
    pca_out_ind <- factoextra::get_pca_ind(pca_out)
    pca_ind_coord <- data.table::as.data.table(pca_out_ind$coord) #coordinates of each point for PCA plot
    
    dat.merge.ind <- cbind(dat, pca_ind_coord)
    data.table::fwrite(dat.merge.ind,
                       "PCA-individuals.csv")
  }
  
  # Create PCA plot with individuals (with groups)
  if (plot.ind.group == TRUE && !is.null(group.ind)) {
    col.factor <- as.factor(dat[, group.ind])
    
    # Set group colours
    if (colour.group == "jet") {
      group.colour.scheme <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
    }
    
    if (colour.group == "spectral") {
      group.colour.scheme <- rev(RColorBrewer::brewer.pal(length(levels(col.factor)),"Spectral"))
    }
    
    if (colour.group == "viridis") {
      group.colour.scheme <- c(viridis::viridis_pal(option = "viridis")(length(levels(col.factor))))
    }
    
    if (colour.group == "inferno") {
      group.colour.scheme <- c(viridis::viridis_pal(option = "inferno")(length(levels(col.factor))))
    }
    
    if (colour.group == "magma") {
      group.colour.scheme <- c(viridis::viridis_pal(option = "magma")(length(levels(col.factor))))
    }
    
    pca_out_ind_group <- factoextra::fviz_pca_ind(pca_out,
                                                  col.ind = col.factor, # color by groups
                                                  palette = group.colour.scheme,
                                                  addEllipses = TRUE, # Concentration ellipses
                                                  ellipse.type = ellipse.type,
                                                  ellipse.level = ellipse.level,
                                                  legend.title = "Groups",
                                                  repel = repel,
                                                  geom = plot.ind.label,
                                                  pointsize = pointsize.group
                                                  )
    
    ggplot2::ggsave(pca_out_ind_group,
                    filename = "PCA plot-ind-groups.pdf",
                    path = path,
                    useDingbats = FALSE)
    
  }
  
  # Create PCA plot with variables
  if (plot.variables == TRUE) {
    # Set colour scheme for variables
    if (colour.var == "solid") {
      var.colour.scheme <- c("#00AFBB", "#E7B800", "#FC4E07")
    }
    
    if (colour.var == "jet") {
      var.colour.scheme <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
    }
    
    if (colour.var == "spectral") {
      var.colour.scheme <- rev(RColorBrewer::brewer.pal(11,"Spectral"))
    }
    
    if (colour.var == "viridis") {
      var.colour.scheme <- c(viridis::viridis_pal(option = "viridis")(50))
    }
    
    if (colour.var == "inferno") {
      var.colour.scheme <- c(viridis::viridis_pal(option = "inferno")(50))
    }
    
    if (colour.var == "magma") {
      var.colour.scheme <- c(viridis::viridis_pal(option = "magma")(50))
    }
    
    if (colour.var == "BuPu") {
      var.colour.scheme <- c(colorRampPalette(RColorBrewer::brewer.pal(9, "BuPu"))(31))
    }
    
    loading_plot <- factoextra::fviz_pca_var(pca_out,
                                             col.var = "contrib",
                                             gradient.cols = var.colour.scheme,
                                             repel = repel,
                                             select.var = list(contrib = var.numb) #plots the top number of contributors
                                             )
    
    # Save loadings plot
    ggplot2::ggsave(loading_plot,
           filename = "PCA plot-variables.pdf",
           path = path,
           useDingbats = FALSE)
    
    # Extract coordinates of variables
    pca_out_var <- factoextra::get_pca_var(pca_out)
    pca_var_coord <- as.data.frame(pca_out_var$coord) #coordinates of each point for PCA plot
    rownames(pca_var_coord) <- colnames(dat[, use.cols])
    
    data.table::fwrite(pca_var_coord,
                       "PCA-variables.csv",
                       row.names = TRUE)
    
  }
  
  # Create combined PCA plot with individuals and variables
  if (plot.combined == TRUE) {
    pca_out_comb <- factoextra::fviz_pca_biplot(pca_out,
                                col.var = "Grey", # Variables color
                                col.ind = "Black",  # Individuals color
                                repel = repel,
                                geom = plot.ind.label,
                                pointsize = pointsize.ind,
                                select.var = list(contrib = var.numb) #plots the top number of contributors
                                )
    
    ggplot2::ggsave(pca_out_comb,
                    filename = "PCA plot-combined.pdf",
                    path = path,
                    useDingbats = FALSE)
  }
  
}



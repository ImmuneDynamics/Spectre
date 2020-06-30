#' Run the PCA algorithm (using stats::prcomp)
#' 
#' Method to run a PCA dimensionality reduction algorithm.
#' A principal component analysis (PCA) is capable of reducing the number of dimensions (i.e. parameters) with minimal effect on the variation of the given dataset.
#' This function will run a PCA calculation (extremely fast) and generate plots (takes time).
#' For individuals (such as samples or patients), a PCA can group them based on their similarities.
#' A PCA is also capable of ranking variables/parameters (such as markers or cell counts) based on their contribution to the variability across a dataset in an extremely fast manner.
#' In cytometry, this can be useful to identify marker(s) that can be used to differentiate between subset(s) of cells.
#' Uses the base R package "stats" for PCA, "factoextra" for PCA, scree and loading plots, "data.table" for saving .csv files, "ggplot2" for saving plots, "gtools" for rearranging data order.
#' More information on PCA plots can be found here http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/.
#'
#' @param dat NO DEFAULT. data.frame.
#' @param use.cols NO DEFAULT. Vector of numbers, reflecting the columns to use for dimensionality reduction (may not want parameters such as "Time" or "Sample").
#' @param scale DEFAULT = TRUE. A logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place.
#' @param scree.plot DEFAULT = TRUE. Option to create scree plots. Note this will require the input of an elbow point during run. Will save generated scree plot.
#' @param component.loading DEFAULT = TRUE. Option to create plots for each component. Requires scree.plot = TRUE.
#' @param variable.contribution DEFAULT = TRUE. Option to create plot showing the contribution of each variable. Horizontal red line represents the average variable contribution if all variables contributed equally. Requires scree.plot = TRUE.
#' @param plot.individuals DEFAULT = TRUE. Option to create PCA plots on individuals (samples/patients).
#' @param plot.ind.label DEFAULT = "point". Option to add text to PCA plots on individuals as an extra identifier. Use c("point", "text") to include both text and point.
#' @param row.names DEFAULT = NULL. Column (as character) that defines individuals. Will be used to place name on plot.individuals.
#' @param plot.ind.group DEFAULT = FALSE. Option to group inidividuals with ellipses. Must specify column that groups individuals with group.ind.
#' @param group.ind DEFAULT = NULL. Column (as character) that defines groups of individuals. Works with plot.ind.group which must be set to TRUE.
#' @param plot.variables DEFAULT = TRUE. Option to create PCA plots on variables (markers/cell counts).
#' @param plot.combined DEFAULT = TRUE. Option to create a combined PCA plot with both individuals and variables.
#' @param repel DEFAULT = FALSE. Option to avoid overlapping text in PCA plots. Can greatly increase plot time if there is a large number of samples.
#' @param path DEFAULT = getwd(). The location to save plots. By default, will save to current working directory. Can be overidden.
#' 
#' @usage run.pca(dat, use.cols, scale = TRUE, scree.plot = TRUE, component.loading = TRUE, variable.contribution = TRUE, plot.individuals = TRUE, plot.ind.label = "point", row.names = NULL, plot.ind.group = FALSE, group.ind = NULL, plot.variables = TRUE, plot.combined = TRUE, repel = FALSE, path = getwd())
#'
#' @examples
#' # Set directory to save files. By default it will save files at get()
#' setwd("/Users/felixmarsh-wakefield/Desktop")
#' 
#' # Run PCA on demonstration dataset
#' Spectre::run.pca(dat = Spectre::demo.start,
#'                 use.cols = c(5:6,8:9,11:13,16:19,21:30,32),
#'                 repel = TRUE
#'                 )
#' 
#' # Compare between groups
#' Spectre::run.pca(dat = Spectre::demo.start,
#'                  use.cols = c(5:6,8:9,11:13,16:19,21:30,32),
#'                  plot.ind.label = c("point", "text"), #individual cells will be labelled as numbers
#'                  plot.ind.group = TRUE,
#'                  group.ind = "Group"
#'                  )
#'         
#' # When prompted, type in "4" and click enter to continue function (this selects the elbow point based off the scree plot)
#' 
#' @author Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#' @export

run.pca <- function(dat,
                 use.cols,
                 scale = TRUE,
                 scree.plot = TRUE,
                 component.loading = TRUE,
                 variable.contribution = TRUE,
                 plot.individuals = TRUE,
                 plot.ind.label = "point",
                 row.names = NULL,
                 plot.ind.group = FALSE,
                 group.ind = NULL,
                 plot.variables = TRUE,
                 plot.combined = TRUE,
                 repel = FALSE,
                 path = getwd()
                 ) {
  
  ## Check that necessary packages are installed
  if(!is.element('stats', installed.packages()[,1])) stop('stats is required but not installed')
  if(!is.element('factoextra', installed.packages()[,1])) stop('factoextra is required but not installed')
  if(!is.element('ggplot2', installed.packages()[,1])) stop('ggplot2 is required but not installed')
  if(!is.element('gtools', installed.packages()[,1])) stop('gtools is required but not installed')
  if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')
  if(!is.element('ggpubr', installed.packages()[,1])) stop('ggpubr is required but not installed')
  
  ## Require packages
  require(stats)
  require(factoextra)
  require(ggplot2)
  require(gtools)
  require(data.table)
  require(ggpubr)
  
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
           path = path)
    
    if (component.loading == TRUE) {
      # Set parameters for plots
      data.loadings <- unclass(pca_out$rotation) #'princomp' equivalent to 'loadings'
      data.load.order <- as.data.frame(data.loadings[gtools::mixedsort(row.names(data.loadings)), ])
      
      # Create bargraphs for PC (up to the elbow)
      max_y_value <- max(data.load.order[,1:elbow.point], na.rm = TRUE)
      max_y_value_p10 <- max_y_value*1.1
      
      min_y_value <- min(data.load.order[,1:elbow.point], na.rm = TRUE)
      min_y_value_p10 <- min_y_value*1.1
      
      # Creates loadings plots for each of the components/dimensions (based on elbow point)
      for (i in 1:elbow.point) {
        p_loading <- ggplot2::ggplot(data = data.load.order, ggplot2::aes(x = rownames(data.load.order), y = data.load.order[,i])) +
          ggplot2::geom_bar(stat="identity", color = "black", fill = "black") +
          #theme_minimal() + #includes grids
          ggplot2::theme_classic() + #removes all grids
          ggplot2::scale_x_discrete(limits = c(rownames(data.load.order))) + #sets order
          ggplot2::scale_y_continuous(limits = c(min_y_value_p10, max_y_value_p10)) +
          ggplot2::labs(title = paste0("PC", i), x = "Variable", y = "Component loading") +
          ggplot2::theme(legend.position = "none", # can be "left" "right" "top" "bottom" "none
                axis.text.x = ggplot2::element_text(colour="black",size=12,angle=90,hjust=1,vjust=1,face="bold"),
                axis.text.y = ggplot2::element_text(colour="black",size=12,angle=0,hjust=1,vjust=0,face="bold"),  
                axis.title.x = ggplot2::element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0,face="bold"),
                axis.title.y = ggplot2::element_text(colour="black",size=12,angle=90,hjust=.5,vjust=1,face="bold"),
                plot.title = ggplot2::element_text(lineheight=.8, face="bold", hjust = 0, size = 18), # hjust = 0.5 to centre
                axis.line = ggplot2::element_line(colour = 'black', size = 0.5),
                axis.ticks = ggplot2::element_line(colour = "black", size = 0.5)
          )
        
        ggplot2::ggsave(p_loading,
                        filename = paste0("PC", i, sep = "_", "loading.pdf"),
                        width = 16,
                        height = 5,
                        path = path)
      }
    }
    
    if (variable.contribution == TRUE) {
      # Creates plot with contributions of each variable
      pca.contrib <- factoextra::fviz_contrib(pca_out, #red dashed line represents the expected average row contributions if the contributions were uniform: 1/nrow(markers)
                                  choice = "var", #variable ("var") or individual ("ind")
                                  axes = 1:elbow.point) #axes refer to PCA number
      
      ggplot2::ggsave(pca.contrib,
             filename = "pca_contribution.pdf",
             width = 16,
             height = 5,
             path = path)
    }
    
  }
  
  # Create PCA plot with individuals
  if (plot.individuals == TRUE) {
    pca_plot_ind <- factoextra::fviz_pca_ind(pca_out,
                                             repel = repel,
                                             geom = plot.ind.label)
    
    ggplot2::ggsave(pca_plot_ind,
                    filename = "PCA plot-individuals.pdf",
                    path = path)
    
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
    
    pca_out_ind_group <- factoextra::fviz_pca_ind(pca_out,
                             col.ind = col.factor, # color by groups
                             # palette = c("#00AFBB",  "#FC4E07"),
                             addEllipses = TRUE, # Concentration ellipses
                             ellipse.type = "confidence",
                             legend.title = "Groups",
                             repel = repel,
                             geom = plot.ind.label
                             )
    
    ggplot2::ggsave(pca_out_ind_group,
                    filename = "PCA plot-ind-groups.pdf",
                    path = path)
    
  }
  
  # Create PCA plot with variables
  if (plot.variables == TRUE) {
    loading_plot <- factoextra::fviz_pca_var(pca_out,
                                      col.var = "contrib",
                                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                      repel = repel)
    
    # Saves loadings plot
    ggplot2::ggsave(loading_plot,
           filename = "Loading plot.pdf",
           path = path)
    
    # Extract coordinates of variables
    pca_out_var <- factoextra::get_pca_var(pca_out)
    pca_var_coord <- as.data.frame(pca_out_var$coord) #coordinates of each point for PCA plot
    rownames(pca_var_coord) <- colnames(dat[, use.cols])
    
    data.table::fwrite(pca_var_coord,
                       "PCA-variables.csv",
                       row.names = TRUE)
    
    data.loadings <- unclass(pca_out$rotation) #shows the loadings values for everything!
    data.table::fwrite(x = as.data.frame(data.loadings),
           file = "loadings.csv",
           row.names = TRUE)
  }
  
  # Create combined PCA plot with individuals and variables
  if (plot.combined == TRUE) {
    pca_out_comb <- factoextra::fviz_pca_biplot(pca_out,
                                col.var = "Grey", # Variables color
                                col.ind = "Black",  # Individuals color
                                repel = repel,
                                geom = plot.ind.label
                                )
    
    ggplot2::ggsave(pca_out_comb,
                    filename = "PCA plot-combined.pdf",
                    path = path)
  }
  
}



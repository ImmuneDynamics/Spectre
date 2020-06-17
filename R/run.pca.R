#' Run the PCA algorithm (using stats::princomp())
#' 
#' Method to run a PCA dimensionality reduction algorithm.
#' A principal component analysis (PCA) is capable of reducing the number of dimensions (i.e. parameters) with minimal effect on the variation of the given dataset.
#' A PCA is also capable of ranking parameters based on their contribution to the variability across a dataset in an extremely fast manner.
#' In cytometry, this can be useful to identify marker(s) that can be used to differentiate between subset(s) of cells.
#' Uses the base R package "stats" for PCA, "factoextra" for PCA, scree and loading plots, "data.table" for saving .csv files, "ggplot2" for saving plots, "gtools" for rearranging data order, "dplyr" for selecting top n values.
#'
#' @param dat NO DEFAULT. data.frame.
#' @param use.cols NO DEFAULT. Vector of numbers, reflecting the columns to use for dimensionality reduction (may not want parameters such as "Time" or "Sample").
#' @param cor DEFAULT = TRUE. A logical value indicating whether the calculation should use the correlation matrix or the covariance matrix. (The correlation matrix can only be used if there are no constant variables.).
#' @param scores DEFAULT = TRUE. A logical value indicating whether the score on each principal component should be calculated.
#' @param pca.plot DEFAULT = TRUE. Option to create PCA plot.
#' @param scree.plot DEFAULT = TRUE. Option to create scree plots. Note this will require the input of an elbow point during run. Will save generated scree plot.
#' @param component.loading DEFAULT = TRUE. Option to create plots for each component. Requires scree.plot = TRUE.
#' @param marker.contribution DEFAULT = TRUE. Option to create plot showing the contribution of each marker. Horizontal red line represents the average marker contribution if all markers contributed equally. Requires scree.plot = TRUE.
#' @param loading.plot DEFAULT = TRUE. Option to create scree plots. Will save generated loading plot.
#' @param individual.samples DEFAULT = FALSE. Option to run above plots on a per sample basis. Only samples that have a cell number greater than the number of parameters/markers will be included.
#' @param sample.code DEFAULT = "FileName". Parameter to define column that will differentiate between samples. Must be quoted. Requires scree.plot = TRUE. Used for "individual.samples".
#' @param top.tally DEFAULT = 10. Number of top markers that contributed to PCA across each sample. Used for "individual.samples".
#' @param path DEFAULT = getwd(). The location to save plots. By default, will save to current working directory. Can be overidden.
#' 
#' @usage run.pca(dat, use.cols, cor = TRUE, scores = TRUE, pca.plot = TRUE, scree.plot = TRUE, component.loading = TRUE, marker.contribution = TRUE, loading.plot = TRUE, individual.samples = FALSE, sample.code = "FileName", top.tally = 10, path = getwd())
#'
#' @examples
#' # Set directory to save files. By default it will save files at get()
#' setwd("/Users/felixmarsh-wakefield/Desktop")
#' 
#' # Run PCA on demonstration dataset
#' Spectre::run.pca(dat = Spectre::demo.start,
#'         use.cols = c(5:6,8:9,11:13,16:19,21:30,32))
#'         
#' # When prompted, type in "4" and click enter to continue function (this selects the elbow point based off the scree plot)
#' 
#' @author Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#' @export

run.pca <- function(dat,
                 use.cols,
                 cor = TRUE,
                 scores = TRUE,
                 pca.plot = TRUE,
                 scree.plot = TRUE,
                 component.loading = TRUE,
                 marker.contribution = TRUE,
                 loading.plot = TRUE,
                 individual.samples = FALSE,
                 sample.code = "FileName",
                 top.tally = 10,
                 path = getwd()
                 ){
  
  ## Check that necessary packages are installed
  if(!is.element('stats', installed.packages()[,1])) stop('stats is required but not installed')
  if(!is.element('factoextra', installed.packages()[,1])) stop('factoextra is required but not installed')
  if(!is.element('ggplot2', installed.packages()[,1])) stop('ggplot2 is required but not installed')
  if(!is.element('gtools', installed.packages()[,1])) stop('gtools is required but not installed')
  if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')
  if(!is.element('dplyr', installed.packages()[,1])) stop('dplyr is required but not installed')
  if(!is.element('ggpubr', installed.packages()[,1])) stop('ggpubr is required but not installed')
  
  ## Require packages
  require(stats)
  require(factoextra)
  require(ggplot2)
  require(gtools)
  require(data.table)
  require(dplyr)
  require(ggpubr)
  
  ## Check dat is a data.table
  if (!data.table::is.data.table(dat)) {
    dat <- data.table::as.data.table(dat)
  }
  
  ## Run PCA
  pca_out <- stats::princomp(as.matrix(dat[, ..use.cols]),
                             cor = cor,
                             scores = scores
                             )
  
  if (pca.plot == TRUE) {
    pca_plot <- factoextra::fviz_pca_ind(pca_out, geom="point")
    
    ggplot2::ggsave(pca_plot,
                    filename = "PCA plot.pdf",
                    path = path)
  }
  
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
      data.loadings <- unclass(pca_out$loadings)
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
          ggplot2::labs(title = paste0("PC", i), x = "Marker", y = "Component loading") +
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
                        width = 12,
                        height = 5,
                        path = path)
      }
    }
    
    if (marker.contribution == TRUE) {
      # Creates plot with contributions of each marker
      pca.contrib <- factoextra::fviz_contrib(pca_out, #red dashed line represents the expected average row contributions if the contributions were uniform: 1/nrow(markers)
                                  choice = "var", #variable ("var") or individual ("ind")
                                  axes = 1:elbow.point) #axes refer to PCA number
      
      ggplot2::ggsave(pca.contrib,
             filename = "pca_contribution.pdf",
             path = path)
    }
    
  }
  
  if (loading.plot == TRUE) {
    loading_plot <- factoextra::fviz_pca_var(pca_out,
                                      col.var = "contrib",
                                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
    
    # Saves loadings plot
    ggplot2::ggsave(loading_plot,
           filename = "Loading plot.pdf",
           path = path)
    
    data.loadings <- unclass(pca_out$loadings) #shows the loadings values for everything!
    data.table::fwrite(x = as.data.frame(data.loadings),
           file = "loadings.csv",
           row.names = TRUE)
  }
  
  if (individual.samples == TRUE) {
    # Create folder
    dir.create(paste0(path, "/samples"), showWarnings = FALSE)
    setwd(paste0(path, "/samples"))
    PrimaryDirectory <- getwd()
    
    # Define sample differentiator
    sample.name <- unique(dat[[sample.code]])
    
    # Create list that will hold contribution data for each sample
    contrib.list <- list()
    
    for (a in sample.name) {
      data.sample <- dat[dat[[sample.code]] == a,]
      
      data.sample.select <- data.sample[,..use.cols]
      as.matrix(colnames(data.sample.select))
      
      if (nrow(data.sample.select) > ncol(data.sample.select)) { #PCA cannot run if there are more markers/columns than cells/rows
        # Run PCA on a single sample
        data.pca.sample <- princomp(data.sample.select, cor = TRUE, scores = TRUE)
        
        setwd(PrimaryDirectory)
        dir.create(paste0(PrimaryDirectory, "/", a), showWarnings = FALSE)
        setwd(paste0(PrimaryDirectory, "/", a))
        OutputDirectory <- getwd()
        
        ## Scree plots
            scree_plot_sample <- factoextra::fviz_eig(data.pca.sample, addlabels = TRUE) #creates scree plot; addlabels adds % to plot
        
            # Saves scree plot
            ggplot2::ggsave(scree_plot_sample,
                            filename = "Scree plot.pdf",
                            path = path)
            
        ## PCA plots
            pca_plot_sample <- factoextra::fviz_pca_ind(data.pca.sample, geom="point")
            
            # Saves PCA plot
            ggplot2::ggsave(pca_plot_sample,
                            filename = "PCA plot.pdf",
                            path = OutputDirectory)
            
        ## Loading plots
            loading_plot_sample <- factoextra::fviz_pca_var(data.pca.sample,
                                                     col.var = "contrib",
                                                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
            
            # Saves loadings plot
            ggplot2::ggsave(loading_plot_sample,
                            filename = "Loading plot.pdf",
                            path = OutputDirectory)
            
            data.loadings.sample <- unclass(data.pca.sample$loadings) #shows the loadings values for everything!
            data.table::fwrite(x = as.data.frame(data.loadings.sample),
                               file = "loadings.csv",
                               row.names = TRUE)
            
        ## Component loading plots
            # Set parameters for plots
            data.loadings.sample <- unclass(data.pca.sample$loadings)
            data.load.order.sample <- as.data.frame(data.loadings.sample[gtools::mixedsort(row.names(data.loadings.sample)), ])
            
            # Create bargraphs for PC (up to the elbow)
            max_y_value <- max(data.load.order.sample[,1:elbow.point], na.rm = TRUE)
            max_y_value_p10 <- max_y_value*1.1
            
            min_y_value <- min(data.load.order.sample[,1:elbow.point], na.rm = TRUE)
            min_y_value_p10 <- min_y_value*1.1
            
            # Creates loadings plots for each of the components/dimensions (based on elbow point)
            for (i in 1:elbow.point) {
              p_loading_sample <- ggplot2::ggplot(data = data.load.order.sample, ggplot2::aes(x = rownames(data.load.order.sample), y = data.load.order.sample[,i])) +
                ggplot2::geom_bar(stat="identity", color = "black", fill = "black") +
                #theme_minimal() + #includes grids
                ggplot2::theme_classic() + #removes all grids
                ggplot2::scale_x_discrete(limits = c(rownames(data.load.order.sample))) + #sets order
                ggplot2::scale_y_continuous(limits = c(min_y_value_p10, max_y_value_p10)) +
                ggplot2::labs(title = paste0("PC", i), x = "Marker", y = "Component loading") +
                ggplot2::theme(legend.position = "none", # can be "left" "right" "top" "bottom" "none
                               axis.text.x = ggplot2::element_text(colour="black",size=12,angle=90,hjust=1,vjust=1,face="bold"),
                               axis.text.y = ggplot2::element_text(colour="black",size=12,angle=0,hjust=1,vjust=0,face="bold"),  
                               axis.title.x = ggplot2::element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0,face="bold"),
                               axis.title.y = ggplot2::element_text(colour="black",size=12,angle=90,hjust=.5,vjust=1,face="bold"),
                               plot.title = ggplot2::element_text(lineheight=.8, face="bold", hjust = 0, size = 18), # hjust = 0.5 to centre
                               axis.line = ggplot2::element_line(colour = 'black', size = 0.5),
                               axis.ticks = ggplot2::element_line(colour = "black", size = 0.5)
                )
              
              ggplot2::ggsave(p_loading_sample,
                              filename = paste0("PC", i, sep = "_", "loading.pdf"),
                              width = 12,
                              height = 5,
                              path = OutputDirectory)
            }
            
        ## Variables (i.e. markers)
            var <- factoextra::get_pca_var(data.pca.sample) #gets data for variables
            
            data.table::fwrite(x = as.data.frame(var$contrib),
                   file = "contributions.csv",
                   row.names = TRUE)
            
            # Visualisation
            pca.contrib.sample <- factoextra::fviz_contrib(data.pca.sample, #red dashed line represents the expected average row contributions if the contributions were uniform: 1/nrow(markers) (1/34 = 2.9%)
                                        choice = "var", #variable ("var") or individual ("ind")
                                        axes = 1:elbow.point) #axes refer to PCA number
            
            ggplot2::ggsave(pca.contrib.sample, filename = "pca_contribution.pdf")
            
        ## Save contribution for later use
            contrib.list[[a]] <- pca.contrib.sample$data
      }
      
      ## Calculate tally for top contributing markers
      setwd(PrimaryDirectory)
      
      contrib.list.top <- lapply(contrib.list, dplyr::top_n, top.tally) #select top n markers that contribute to variation across select PCA/dimensions
      contrib.freq <- as.data.frame(table(data.table::rbindlist(contrib.list.top)$name)) #calculates frequency of markers that were in the top n contributors of variation
      
      p.contrib <- ggpubr::ggbarplot(contrib.freq, x = "Var1", y = "Freq", fill = "steelblue", 
                                     color = "steelblue", sort.val = "desc", top = Inf, main = "Occurrences in top 10 contributions", 
                                     xlab = "Marker", ylab = paste0("Occurrence (out of ", length(contrib.list.top), ")"), xtickslab.rt = 45, 
                                     ggtheme = ggplot2::theme_minimal(), sort.by.groups = FALSE) +
        ggplot2::geom_text(ggplot2::aes(label = Freq), vjust = -0.5) #adds values on top of bar
      
      ggplot2::ggsave(p.contrib,
                      filename = paste0("Top ", top.tally, " contributors tally.pdf"),
                      path = PrimaryDirectory)
    }
  }
  
}



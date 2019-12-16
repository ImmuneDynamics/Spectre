#' run.pca - ...
#'
#' @usage run.pca(x, use.cols, ...)
#'
#' @param x data.frame. No default.
#' @param use.cols Vector of numbers, reflecting the columns to use for clustering. No default.
#' @param cor A logical value indicating whether the calculation should use the correlation matrix or the covariance matrix. (The correlation matrix can only be used if there are no constant variables.). Default = TRUE.
#' @param scores A logical value indicating whether the score on each principal component should be calculated. Default = TRUE.
#' @param scree.plot Option to create scree plots. Note this will require the input of an elbow point during run. Will save generated scree plot. Default = TRUE.
#' @param component.loading Option to create plots for each component. Requires scree.plot = TRUE. Default = TRUE.
#' @param marker.contribution Option to create plot showing the contribution of each marker. Horizontal red line represents the average marker contribution if all markers contributed equally. Requires scree.plot = TRUE. Default = TRUE.
#' @param loading.plot Option to create scree plots. Will save generated loading plot. Default = TRUE.
#' @param individual.samples Option to run above plots on a per sample basis. Only samples that have a cell number greater than the number of parameters/markers will be included. Default = FALSE.
#' @param sample.code Parameter to define column that will differentiate between samples. Must be quoted. Requires scree.plot = TRUE.  Used for "individual.samples".Default = "FileName"
#' @param top.tally Number of top markers that contributed to PCA across each sample. Used for "individual.samples". Defalut = 10.
#' @param path Location to save plots. Default = getwd() (working directory)
#'
#'This function runs a principal component analysis (PCA) on a dataframe with cells (rows) vs markers (columns), returning chosen figures. Uses the base R package "stats" for PCA, "factoextra" for scree and loading plots, "data.table" for saving .csv files, "ggplot2" for saving plots, "gtools" for rearranging data order, "dplyr" for selecting top n values.
#'
#' @export

run.pca <- function(x,
                 use.cols,
                 cor = TRUE,
                 scores = TRUE,
                 scree.plot = TRUE,
                 component.loading = TRUE,
                 marker.contribution = TRUE,
                 loading.plot = TRUE,
                 individual.samples = FALSE,
                 sample.code = "FileName",
                 top.tally = 10,
                 path = getwd()
                 ){
  
  pca_out <- stats::princomp(as.matrix(x[use.cols]),
                             cor = cor,
                             scores = scores,
                             scree.plot = scree.plot,
                             component.loading = component.loading,
                             marker.contribution = marker.contribution,
                             loading.plot = loading.plot,
                             individual.samples = individual.samples,
                             sample.code = sample.code,
                             top.tally = top.tally,
                             path = path
                             )
  
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
    sample.name <- unique(x[[sample.code]])
    
    # Create list that will hold contribution data for each sample
    contrib.list <- list()
    
    for (a in sample.name) {
      data.sample <- x[x[[sample.code]] == a,]
      
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



#' create.dt - convert a FlowFrame or Seurat or SingleCellExperiment to a data.table
#'
#' This function converts a Seurat or SingleCellExperiment object into a list containing a data.table, with vectors of gene and dimensionality reduction
#'
#' @param dat NO DEFAULT. A Seurat or SingleCellExperiment object.
#' @param from DEFAULT = NULL. By default, the class of object will be detected automatically, but this can be overwritten using from. Can be from = 'Seurat' or 'SingleCellExperiment'.
#' 
#' @usage create.dt(dat, from)
#'
#' @author
#' Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://github.com/ImmuneDynamics/Spectre}.
#' 
#' @import data.table
#' 
#' @export

create.dt <- function(dat,
                      from = NULL)
{
  
  ######################################################################################################
  ### Setup
  ######################################################################################################
      
  ### Packages
  
      require('Spectre')
      require('data.table')
  
  ### Determine class of object
      
      if(class(dat)[1] == 'Seurat'){
        object.type <- "Seurat"
      }
      
      if(class(dat)[1] == 'SingleCellExperiment'){
        object.type <- "SingleCellExperiment"
      }
      
      if(class(dat)[1] == 'flowFrame'){
        object.type <- "flowFrame"
      }
      
  ### Optional overview Overwrite
      
      if(!is.null(from)){
        object.type <- from
      }
      
  ### Warning
      
      if(!exists('object.type')){
        stop("Could not determine the object type. Currently only Seurat objects are supported. You can also try manually specifying the object type using the 'from' argument (e.g. from = 'Seurat'")
      }
  
  ######################################################################################################
  ### OPTION: Seurat objects
  ######################################################################################################
  
      if(object.type == "Seurat"){
        
        message(object.type, ' detected')
        
        ### Packages
        
            require('dplyr')
            require('Seurat')
            require('patchwork')
        
        ### Extract
        
            a <- GetAssayData(object = dat)
            assays <- names(dat@assays)
            dim.reds <- names(dat@reductions)

            var.features <- VariableFeatures(dat)
            
            if(length(var.features) > 0){
              var.features.top10 <- head(VariableFeatures(dat), 10)
            }
            

        ### Genes and cells
            
            geneNames <- a@Dimnames[[1]]
            cellNames <- a@Dimnames[[2]]

        ### Start data.table
            
            res <- as.data.table(cellNames)
            names(res) <- 'cellNames'
            
        ### Add data from slots
   
            for(i in assays){
              # i <- assays[[1]]

              types <- vector()
              
              if(ncol(dat@assays[[i]]@counts) > 0){
                types <- c(types, 'counts')
                
                x1 <- GetAssayData(object = dat, assay = i, slot = 'counts')
                x1 <- as.data.table(x1)
                x1 <- data.table::transpose(x1)
                names(x1) <- geneNames 
                names(x1) <- paste0(names(x1), "_", i, "_", 'counts')
                res <- cbind(res, x1)
              }
              
              if(ncol(dat@assays[[i]]@data) > 0){
                types <- c(types, 'data')
                
                x2 <- GetAssayData(object = dat, assay = i, slot = 'data')
                x2 <- as.data.table(x2)
                x2 <- data.table::transpose(x2)
                names(x2) <- geneNames 
                names(x2) <- paste0(names(x2), "_", i, "_", 'data')
                res <- cbind(res, x2)
              }
              
              if(ncol(dat@assays[[i]]@scale.data) > 0){
                types <- c(types, 'scale.data')
                
                x3 <- GetAssayData(object = dat, assay = i, slot = 'scale.data')
                x3 <- as.data.table(x3)
                x3 <- data.table::transpose(x3)
                names(x3) <- geneNames 
                names(x3) <- paste0(names(x3), "_", i, "_", 'scale.data')
                res <- cbind(res, x3)
              }
              
              rm(i)
              rm(x1)
              rm(x2)
              rm(x3)
            }
            
        ### Add dim reductions
            
            for(i in dim.reds){
              # i <- dim.reds[[2]]
              
              tmp <- dat@reductions[[i]]@cell.embeddings
              tmp <- as.data.table(tmp)
              
              names(tmp) <- paste0(i, '_', names(tmp))
              res <- cbind(res, tmp)
            }
            
        ### Wrap up and return
 
            final.res <- list()
            
            final.res$data.table <- res
            final.res$geneNames <- geneNames
            final.res$cellNames <- cellNames
            
            if(length(var.features) > 0){
              final.res$var.features <- var.features
              final.res$var.features.top10 <- var.features.top10
            }
            
            final.res$assays <- paste0('_', assays)
            final.res$slots <- paste0('_', types)
            final.res$dim.reds <- paste0(dim.reds, '_')

            message(paste0("Converted a ", object.type, " object into a data.table stored in a list"))
            return(final.res)
      }
      
  ######################################################################################################
  ### OPTION: SingleCellExperiment
  ######################################################################################################
  
      if(object.type == "SingleCellExperiment"){
        
        message(object.type, ' detected')
        
        ### Packages
        
            require('SingleCellExperiment')
        
        ###
        
            geneNames <- rownames(dat) # genes
            geneNames
            
            cellNames <- colnames(dat) # cells
            cellNames
        
        ###
            assays <- names(dat@assays@data)
            assays
            
            dim.reds <- names(reducedDims(dat))
            dim.reds
        
        ###
            res <- as.data.table(cellNames)
            names(res) <- 'cellNames'
        
        ### Add metadata
            
            if(exists('dat@colData')){
              col.meta <- dat@colData
              col.meta <- as.data.table(col.meta)
              meta.cols <- names(col.meta)
              
              res <- cbind(res, col.meta)
              res
            }

        ### Assay data
            
            for(i in assays){
              # i <- assays[[1]]
              
              tmp <- dat@assays@data[[i]]
              # dim(tmp)
              # colnames(tmp) # columns = cells
              tmp <- as.data.table(tmp)
              tmp <- data.table::transpose(tmp)  
              names(tmp) <- geneNames
              names(tmp) <- paste0(names(tmp), "_", i)
              
              res <- cbind(res, tmp)
              
              rm(i)
              rm(tmp)
            }
        
        ### DimRed data
            
            for(i in dim.reds){
              # i <- dim.reds[[2]]
              
              tmp <- reducedDims(dat)[[i]]
              # dim(tmp)
              # colnames(tmp) # columns = DRs
              tmp <- as.data.table(tmp) # Dont transpose
              
              new.names <- paste0(i, "_", c(1:length(names(tmp))))
              names(tmp) <- new.names
              
              res <- cbind(res, tmp)
              
              rm(i)
              rm(tmp)
            }
            
        ### Wrap up and return
            
            final.res <- list()
            
            final.res$data.table <- res
            final.res$geneNames <- geneNames
            final.res$cellNames <- cellNames
            final.res$assays <- paste0('_', assays)
            final.res$dim.reds <- paste0(dim.reds, '_')
            
            message(paste0("Converted a ", object.type, " object into a data.table stored in a list"))
            return(final.res)
      }

}


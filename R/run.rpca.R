#' Run rPCA for batch correction
#' 
#' Use Seurat rPCA to do batch correction.
#' For more details on rPCA, see https://satijalab.org/seurat/articles/integration_rpca.html.
#'
#' @param dat NO DEFAULT. A data.table
#' @param use.cols NO DEFAULT. 
#' A vector of character column names to apply batch correction to.
#' @param batch.col Character. The column in the data in data.table that identifies
#' which batch each cell belongs to.
#' @param reference DEFAULT NULL. Whether to align to batches to a given batch.
#' If yes, then supply this parameter with the name of the batch you want to align the other batches to.
#' @param k.anchor DEFAULT 5. Passed to Seurat's FindIntegrationAnchors function.
#' Essentially, it determines the number of neighbors (k) to use when 
#' Seurat's FindIntegrationAnchors is picking anchors.
#'
#' @return batch corrected data.table
#'
#' 
#' @export
#' 

run.rpca <- function(dat,
                     use.cols,
                     batch.col,
                     reference = NULL,
                     k.anchor = 5
){
  
  ### Setup
  
      require('data.table')
      require('Matrix')
      require('Seurat')
      require('SeuratObject')
      
  ### Test data
  
      # use.cols <- cluster.cols
      # batch.col <- 'Batch'
      # dat <- cell.dat
      # reference = NULL
      # k.anchor = 5
  
  ### Split data by batch/dataset and create Seurat objects
      
      message('Running rPCA')
      message(' -- setting up data (1/5)')
      
      dat.start <- dat
      raw.cols <- use.cols
      
      new.names <- paste0('Col', c(1:length(raw.cols)))
      match <- data.table('Start' = raw.cols, 'New' = new.names)
      dat <- dat[,c(batch.col, raw.cols), with = FALSE]
      names(dat)[c(2:length(names(dat)))] <- new.names
      dat$CELLBARCODE <- c(1:nrow(dat))
      batches <- unique(dat[[batch.col]])
      
      org.list <- list()
      res.list <- list()
      
      for(i in batches){
        # i <- batches[[1]]
        
        message(paste0('    ...', i))
        
        rws <- dat[[batch.col]] == i
        
        rw.org <- dat[rws, 'CELLBARCODE', with = FALSE]
        org.list[[i]] <- rw.org
        
        dat.batch <- dat[rws,..new.names]
        
        counts <- as.matrix(dat.batch)
        rownames(counts) <- paste0('Cell-', i, '-', c(1:nrow(counts)))
        counts <- t(counts)
        counts <- Matrix::Matrix(counts, sparse = TRUE)
        
        srt <- CreateSeuratObject(counts = counts, assay = 'cyto')
        VariableFeatures(srt) <- rownames(counts)
        res.list[[i]] <- srt
        
        rm(rws)
      }

      rm(dat)
      gc()
  
  ### Select integration features, scale data, and run PCA
  
      message(' -- performing scaling and PCA (2/5)')
      
      features <- SelectIntegrationFeatures(object.list = res.list)

      for(i in batches){
        
        message(paste0('    ...', i))
        
        res.list[[i]]@assays$cyto@layers$data <- as.matrix(res.list[[i]]@assays$cyto@layers$counts)
        res.list[[i]]@assays$cyto@cells[['data']] <- res.list[[i]]@assays$cyto@cells[['counts']]
        res.list[[i]]@assays$cyto@features[['data']] <- res.list[[i]]@assays$cyto@features[['counts']]
        res.list[[i]] <- Seurat::ScaleData(res.list[[i]], features = features, verbose = FALSE, assay = 'cyto')
        res.list[[i]] <- Seurat::RunPCA(res.list[[i]], features = features, verbose = FALSE, assay = 'cyto')
      }
  
  ### Find integration anchors across datasets
  
      message(' -- finding integration anchors (3/5)')
      
      if(is.null(reference)){
        immune.anchors <- FindIntegrationAnchors(object.list = res.list, 
                                                 anchor.features = features, 
                                                 dims = 1:(length(features)-1), 
                                                 k.anchor = k.anchor,
                                                 reduction = 'rpca')
      }
      
      if(!is.null(reference)){
        immune.anchors <- FindIntegrationAnchors(object.list = res.list, 
                                                 anchor.features = features, 
                                                 dims = 1:(length(features)-1), 
                                                 k.anchor = k.anchor,
                                                 reduction = 'rpca', 
                                                 reference = which(names(res.list) == reference))
      }
  
  ### Integrate the data
  
      message(' -- integrating data (4/5)')
      
      immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:(length(use.cols)-1))
      DefaultAssay(immune.combined) <- "integrated"
  
  ### Re-construct Spectre object
  
      message(' -- constructing final data (5/5)')
      
      ordr <- rbindlist(org.list)
      
      final <- as.data.table(t(as.matrix(immune.combined@assays$integrated@data)))
      final$CELLBARCODE <- ordr
      
      setorderv(final, 'CELLBARCODE')
      final[['CELLBARCODE']] <- NULL
      
      final <- final[,c(match[[2]]), with = FALSE]
      names(final) <- paste0(match[[1]], '_rPCA_aligned')

      dat.start <- cbind(dat.start, final)
  
  ### Return
  
      return(dat.start)
}

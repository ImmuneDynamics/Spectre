#' run rpca
#'
#' @import data.table Matrix Seurat
#'
#' @export

run.rpca <- function(dat,
                     use.cols,
                     batch.col,
                     reference = NULL,
                     k.anchor = 5
){
  
  ### Setup
      
      require('data.table')
      require('Seurat')
      require('Matrix')
  
  ### Demo data
      
      # use.cols <- cluster.cols
      # #type <- 'asinh'
      # batch.col <- 'Batch'
      # 
      # dat <- cell.dat
      
      # dat <- create.spectre(cell.dat, meta.cols = c('Sample', 'Group', 'Batch'))
      # 
      # as.matrix(names(dat@data$raw))
      # 
      # dat@data$asinh <- dat@data$raw[,names(dat@data$raw)[c(11:18)], with = FALSE] 
      # dat@data$raw <- dat@data$raw[,names(dat@data$raw)[c(1:8)], with = FALSE]
      # 
      # names(dat@data$asinh) <- gsub('_asinh', '', names(dat@data$asinh))
      # 
      # dat
      # 
  
  ### Checks
      
      # if(class(dat) != 'spectre'){
      #   stop("'dat' must be a 'spectre' object")
      # }
  
  ### Setup
      
      message('Running rPCA')
      message(' -- setting up data')
      
      batches <- unique(dat[[batch.col]])
      batches
      
      org.list <- list()
      res.list <- list()
      
  ### Adjust names and subset
      
      clnm <- c(1:length(use.cols))
      clnm <- paste0('Col', sprintf("%05d",clnm))
      
      x <- dat[,use.cols, with = FALSE]
      names(x) <- clnm
      x <- cbind(x, dat[,..batch.col])
      x$CELLBARCODE <- c(1:nrow(x))
      
      start.cols <- use.cols
      use.cols <- clnm
      
      gc()
      
  ### Split data by batch/dataset and create Seurat objects
      
      for(i in batches){
        # i <- batches[[1]]
        
        message(paste0('    ...', i))
        
        rw.org <- x[x[[batch.col]] == i, 'CELLBARCODE', with = FALSE]
        org.list[[i]] <- rw.org
        
        rws <- x[[batch.col]] == i
        
        dat.batch <- x[rws,use.cols, with = FALSE]
        
        counts <- as.matrix(dat.batch)
        
        counts <- Matrix::Matrix(counts, sparse = TRUE)
        rownames(counts) <- paste0('Cell-', i, '-', c(1:nrow(counts)))
        counts <- t(counts)
        
        srt <- CreateSeuratObject(counts = counts, assay = 'cyto')
        srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = length(use.cols))
        res.list[[i]] <- srt
      }
      
      rm(x)
      gc()
      
      use.cols
      new.names <- rownames(res.list[[1]]@assays$cyto@counts)
      
  ### Select integration features, scale data, and run PCA
      
      message(' -- performing scaling and PCA')
      
      features <- SelectIntegrationFeatures(object.list = res.list)
      features
      
      for(i in batches){
        
        message(paste0('    ...', i))
        
        res.list[[i]] <- ScaleData(res.list[[i]], features = features, verbose = FALSE, assay = 'cyto')
        res.list[[i]] <- RunPCA(res.list[[i]], features = features, verbose = FALSE, assay = 'cyto')
        
      }
      
  ### Find integration anchors across datasets
      
      message(' -- finding integration anchors')
      
      
      if(is.null(reference)){
        immune.anchors <- FindIntegrationAnchors(object.list = res.list, 
                                                 anchor.features = features, 
                                                 dims = 1:(length(use.cols)-1), 
                                                 k.anchor = k.anchor,
                                                 reduction = 'rpca')
        
        immune.anchors
        
      }
      
      if(!is.null(reference)){
        immune.anchors <- FindIntegrationAnchors(object.list = res.list, 
                                                 anchor.features = features, 
                                                 dims = 1:(length(use.cols)-1), 
                                                 k.anchor = k.anchor,
                                                 reduction = 'rpca', 
                                                 reference = which(names(res.list) == reference))
        
        immune.anchors
      }
      
      rm(res.list)
      gc()
      
  ### Integrate the data
      
      message(' -- integrating data')
      
      immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:(length(use.cols)-1))
      DefaultAssay(immune.combined) <- "integrated"
  
  ### Re-construct Spectre object
      
      message(' -- constructing final data')
      
      ordr <- rbindlist(org.list)
      rm(org.list)
      
      final <- as.data.table(t(immune.combined@assays$integrated@data))
      rm(immune.combined)
      gc()
      
      final$CELLBARCODE <- ordr
      
      setorderv(final, 'CELLBARCODE')
      final[['CELLBARCODE']] <- NULL
      
  ### Correct names
      
      final <- final[,sort(names(final)), with = FALSE]
      names(final) <- start.cols
      names(final) <- paste0(names(final), '_rPCA_aligned')
      
      dat <- cbind(dat, final)
      
  ### Return
  
      return(dat)
}



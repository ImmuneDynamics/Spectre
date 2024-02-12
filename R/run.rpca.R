#' run.rpca
#'
#' @import data.table Matrix Seurat
#'
#' @export

run.rpca <- function(dat,
                     use.cols,
                     #type, # data.table from the 'data' slot -- 'asinh'
                     batch.col, # Column from the 'meta' table -- Batch,
                     reference = NULL,
                     k.anchor = 5
){
  
  ### Setup
  
  # require: Seurat, data.table, Matrix
  check_packages_installed(c("Seurat", "Matrix"))
  
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
  
  ### Split data by batch/dataset and create Seurat objects
      
      message('Running rPCA')
      message(' -- setting up data (1/5)')
      
      dat.start <- dat
      
      raw.cols <- use.cols

      new.names <- paste0('Col', c(1:length(raw.cols)))
      match <- data.table('Start' = raw.cols, 'New' = new.names)
      match
      
      dat <- dat[,c(batch.col, raw.cols), with = FALSE]
      names(dat)[c(2:length(names(dat)))] <- new.names
      dat
      
      dat$CELLBARCODE <- c(1:nrow(dat))
      
      batches <- unique(dat[[batch.col]])
      batches
      
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
        counts <- Matrix::Matrix(counts, sparse = TRUE)
        rownames(counts) <- paste0('Cell-', i, '-', c(1:nrow(counts)))
        counts <- t(counts)
        
        srt <- CreateSeuratObject(counts = counts, assay = 'cyto')
        srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = length(use.cols))
        res.list[[i]] <- srt
        
        rm(rws)
      }
      
      org.list
      res.list
      
      rm(dat)
      gc()
  
  ### Select integration features, scale data, and run PCA
      
      message(' -- performing scaling and PCA (2/5)')
      
      features <- SelectIntegrationFeatures(object.list = res.list)
      features
      
      for(i in batches){
        
        message(paste0('    ...', i))
        
        res.list[[i]] <- ScaleData(res.list[[i]], features = features, verbose = FALSE, assay = 'cyto')
        res.list[[i]] <- RunPCA(res.list[[i]], features = features, verbose = FALSE, assay = 'cyto')
        
      }
  
  ### Find integration anchors across datasets
      
      message(' -- finding integration anchors (3/5)')
      
      
      if(is.null(reference)){
        immune.anchors <- FindIntegrationAnchors(object.list = res.list, 
                                                 anchor.features = features, 
                                                 dims = 1:(length(features)-1), 
                                                 k.anchor = k.anchor,
                                                 reduction = 'rpca')
        
        immune.anchors
        
      }
      
      if(!is.null(reference)){
        immune.anchors <- FindIntegrationAnchors(object.list = res.list, 
                                                 anchor.features = features, 
                                                 dims = 1:(length(features)-1), 
                                                 k.anchor = k.anchor,
                                                 reduction = 'rpca', 
                                                 reference = which(names(res.list) == reference))
        
        immune.anchors
      }
  
  ### Integrate the data
  
      message(' -- integrating data (4/5)')
      
      immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:(length(use.cols)-1))
      DefaultAssay(immune.combined) <- "integrated"
  
  ### Re-construct Spectre object
      
      message(' -- constructing final data (5/5)')
      
      ordr <- rbindlist(org.list)
      
      final <- as.data.table(t(immune.combined@assays$integrated@data))
      final$CELLBARCODE <- ordr
      
      setorderv(final, 'CELLBARCODE')
      final[['CELLBARCODE']] <- NULL
      
      final <- final[,c(match[[2]]), with = FALSE]
      names(final) <- paste0(match[[1]], '_rPCA_aligned')
      final

      dat.start <- cbind(dat.start, final)
  
  ### Return
  
  return(dat.start)
  
}



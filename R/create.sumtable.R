#' create.sumtables - create a data.table 'summarising' cellular data by sample and population/cluster.
#'
#' This function summarises cellular data and generates a summary data.table
#'
#' @param dat NO DEFAULT. A data.table containing cells (rows) vs features/markers (columns). One column must represent sample names, and another must represent populations/clusters.
#' @param sample.col NO DEFAULT. Character. Name of the sample column (e.g. "Sample").
#' @param pop.col NO DEFAULT. Character. Name of the population/cluster column (e.g. "Population", "Cluster").
#' @param use.cols NO DEFAULT. A character vector indicating the columns to be measured (e.g. cellular columns -- c("CD45", "CD3e") etc).
#' @param annot.cols DEFAULT = NULL. A character vector indicating the columns to be included as annotation columns (e.g. c("Batch", "Group") etc). 
#' @param parent.col DEFAULT = NULL. A character entry indicating a column that represents the 'lineage' each population belongs to (e.g. 'CD4 T cells' may belong to the 'T cells' lineage). Use this to also calculate each population as a percentage of lineage.
#' @param counts DEFAULT = NULL. If you wish to calculate the actual number of cells per sample, a data.table containing the sample names (in column 1) and cell counts per sample (column 2). 
#' @param perc.pos DEFAULT = NULL. If you wish to calculate the percentage of each population that is 'positive' for a marker, you can provide a data.table containing the mark names (in column 1) and cut off values for positivity (column 2). 
#' @param func DEFAULT = "median". Can be "median" or "mean". Defines the type of function for calculating MFI data.
#' @param sep DEFAULT = " -- ". Character separation of the measurement type and the population (e.g. MFI of CD4 -- T cells)
#' 
#' @usage create.sumtable(dat, sample.col, pop.col, use.cols, annot.cols, counts, func, sep)
#'
#' @examples
#' ## Calculate and export results from demonstration data
#' dat <- Spectre::demo.clustered
#' counts <- data.frame('Sample' = unique(dat[['Sample']]), 'Counts' = c(rep(100000, 6), rep(1000000, 6)))
#' 
#' sum.dat <- create.sumtable(dat = dat,
#'                             sample.col = "Sample",
#'                             pop.col = "Population", 
#'                             use.cols = names(dat)[c(11:19)],
#'                             counts = counts)
#' 
#' @author
#' Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#' 
#' @export

create.sumtable <- function(dat,
                            sample.col,
                            pop.col,
                            use.cols,
                            annot.cols = NULL,
                            parent.col = NULL,
                            counts = NULL,
                            perc.pos = NULL,
                            func = 'median',
                            sep = " -- "){
  
  ### Packages
      if(!is.element('Spectre', installed.packages()[,1])) stop('pheatmap is required but not installed')
      if(!is.element('data.table', installed.packages()[,1])) stop('pheatmap is required but not installed')
      
      require(Spectre)
      require(data.table)
  
  ### Demo data
  
      # dat <- Spectre::demo.clustered
      # sample.col <- 'Sample'
      # pop.col <- 'Population'
      # annot.cols <- c('Group', 'Batch')
      # use.cols <- names(dat)[c(11:19)]
      # sep = " -- "
      # func = 'median'
      # 
      # parent.col = 'Lineage'
      # 
      # counts <- data.table('Sample' = unique(dat[['Sample']]),
      #                      'Counts' = c(rep(100000, 6), rep(1000000, 6)))
      # perc.pos <- data.table('Marker' = c('Ly6C_asinh', 'CD11b_asinh'),
      #                        'Cutoff' = c(3, 3.5))
      # 
      # lin.tb <- data.table('Population' = unique(dat$Population),
      #                      'Lineage' = c('Resident' , 'Infiltrating', 'Infiltrating','Infiltrating','Infiltrating','Infiltrating'))
      # 
      # dat <- do.add.cols(dat, 'Population', lin.tb, 'Population')
  
  ### Checks
  
      dat <- as.data.table(dat)
      
      if(!is.null(counts)){
        counts <- as.data.table(counts)
      }
  
  ### Per sample loop
      
      message('Creating summary table')
      
      samps <- sort(unique(dat[[sample.col]]))
      pops <- sort(unique(dat[[pop.col]]))
  
      if(!is.null(parent.col)){
        parents <- sort(unique(dat[[parent.col]]))
      }
      
      res.list <- list()
      
      for(i in samps){
        # i <- samps[[7]]

        message(paste0(' -- processing sample ', i))
        
        ## Initialise table
            dt <- as.data.table(pops)
            names(dt) <- pop.col
        
        ## Subset sample
            temp <- dat[dat[[sample.col]] == i,]
            
            if(is.null(annot.cols)){
              annots <- temp[1,c(sample.col),with = FALSE]  
            }
            
            if(!is.null(annot.cols)){
              annots <- temp[1,c(sample.col, annot.cols),with = FALSE]  
            }
            
        ## Population percentages
            percent <- temp[, .(Percent = .N), by = pop.col]
            percent[[2]] <- percent[[2]] / sum(percent[[2]]) * 100
            names(percent) <- c(pop.col, 'Percent of sample')
            
        ## Cell counts
            if(!is.null(counts)){
              ttl <- counts[counts[[1]] == i,2]
              ttl <- ttl[[1]]
              
              counts.per.sample <- percent
              counts.per.sample[[2]] <- (counts.per.sample[[2]] * ttl) / 100
              names(counts.per.sample) <- c(pop.col, 'Cells per sample')
            }
        
        ## Percent of parent
            
            if(!is.null(parent.col)){
              
              all.perc.of.parent <- list()
            
              for(a in parents){
                # a <- parents[[1]]
                
                tp <- temp[temp[[parent.col]] == a,]
                
                percent.of.parent <- tp[, .(Percent = .N), by = pop.col]
                percent.of.parent[[2]] <- percent.of.parent[[2]] / sum(percent.of.parent[[2]]) * 100
                names(percent.of.parent) <- c(pop.col, paste0('Percent of ', parent.col))
                
                all.perc.of.parent[[a]] <- percent.of.parent
                
                rm(a)
                rm(tp)
                rm(percent.of.parent)
              }
              
              all.perc.of.parent <- rbindlist(all.perc.of.parent, fill = TRUE)
            }
            
        ## MFIs
            
            mfis <- do.aggregate(temp, use.cols = use.cols, by = pop.col, func = func)
            names(mfis)[c(2:length(names(mfis)))] <- paste0('MFI of ', names(mfis)[c(2:length(names(mfis)))])
            
        ## Percent positive

            if(!is.null(perc.pos)){
              
              all.pos.list <- list()
              
              for(o in perc.pos[,1][[1]]){
                # o <- perc.pos[,1][[1]][[1]]
                
                ctf <- perc.pos[perc.pos[[1]] == o,2]
                ctf <- ctf[[1]]
                
                pos.res.lst <- list()
                for(f in pops){
                  # f <- pops[[1]]
                  pop.dt <- temp[temp[[pop.col]] == f,] 
                  
                  if(nrow(pop.dt) != 0){
                    pos.res.lst[[f]] <- nrow(pop.dt[pop.dt[[o]] > ctf,]) / nrow(pop.dt) * 100
                  }
                  
                  if(nrow(pop.dt) == 0){
                    pos.res.lst[[f]] <- NA
                  }
                }
                
                all.pos.list[[o]] <- as.data.table(unlist(pos.res.lst))
              }
              
              perc.res <- data.table('Population' = pops)
              names(perc.res) <- pop.col
              
              perc.res <- cbind(perc.res, as.data.table(all.pos.list))
              
              names(perc.res)[c(2:length(names(perc.res)))] <- paste0("Percent expressing ", names(perc.res)[c(2:length(names(perc.res)))])
            }

        ## Combine results
            
            ## Percent 
            dt <- do.add.cols(dt, pop.col, percent, pop.col, show.status = FALSE)
            
            ## Counts
            if(!is.null(counts)){
              dt <- do.add.cols(dt, pop.col, counts.per.sample, pop.col, show.status = FALSE)
            }
            
            ## Percent of parent
            if(!is.null(parent.col)){
              dt <- do.add.cols(dt, pop.col, all.perc.of.parent, pop.col, show.status = FALSE)
            }

            ## MFIs
            dt <- do.add.cols(dt, pop.col, mfis, pop.col, show.status = FALSE)
            
            ## Percent positive
            if(!is.null(perc.pos)){
              dt <- do.add.cols(dt, pop.col, perc.res, pop.col, show.status = FALSE)
            }
            
        ## Reshape
            
            dt <- melt(dt, id.vars=c(pop.col))
            dt$Measurement <- paste0(dt$variable, sep, dt[[pop.col]])
            dt <- dt[,c('Measurement', 'value'), with = FALSE]
            names(dt) <- c("Measurement", 'Value')

            Tdt <- transpose(dt)
            
            names(Tdt) <- dt[[1]]
            Tdt <- Tdt[-1,]
            
            Tdt[] <- lapply(Tdt, function(x) as.numeric(x))
            Tdt <- as.data.table(cbind(annots, Tdt))
 
            res.list[[i]] <- Tdt
            
            rm(i)
            rm(dt)
            rm(Tdt)
      }
      
  ### Finalise and return
      
      res <- rbindlist(res.list, fill = TRUE) 
      return(res)
      
} 

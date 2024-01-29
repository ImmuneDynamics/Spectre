#' create.sumtable - create a data.table 'summarising' cellular data by sample and population/cluster.
#'
#' This function summarises cellular data and generates a summary data.table
#'
#' @param dat NO DEFAULT. A data.table containing cells (rows) vs features/markers (columns). One column must represent sample names, and another must represent populations/clusters.
#' @param sample.col NO DEFAULT. Character. Name of the sample column (e.g. "Sample").
#' @param pop.col NO DEFAULT. Character. Name of the population/cluster column (e.g. "Population", "Cluster").
#' @param use.cols DEFAULT = NULL A character vector indicating the columns to be measured (e.g. cellular columns -- c("CD45", "CD3e") etc).
#' @param annot.cols DEFAULT = NULL. A character vector indicating the columns to be included as annotation columns (e.g. c("Batch", "Group") etc). 
#' @param parent.col DEFAULT = NULL. A character entry indicating a column that represents the 'lineage' each population belongs to (e.g. 'CD4 T cells' may belong to the 'T cells' lineage). Use this to also calculate each population as a percentage of lineage.
#' @param counts DEFAULT = NULL. If you wish to calculate the actual number of cells per sample, a data.frame or data.table containing the sample names (in column 1) and cell counts per sample (column 2). 
#' @param perc.pos DEFAULT = NULL. If you wish to calculate the percentage of each population that is 'positive' for a marker, you can provide a data.table containing the mark names (in column 1) and cut off values for positivity (column 2). 
#' @param double.pos DEDAULT = NULL. List of vectors, each vector containing the names of multiple markers you wish to calculate % positive for (e.g. CD38+HLADR+). Generates 'and' and 'or' combinations.
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
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://github.com/ImmuneDynamics/Spectre}.
#' 
#' @export

create.sumtable <- function(dat,
                            sample.col,
                            pop.col,
                            use.cols = NULL,
                            annot.cols = NULL,
                            parent.col = NULL,
                            counts = NULL,
                            perc.pos = NULL,
                            double.pos = NULL,
                            func = 'median',
                            sep = " -- "){
  
      # require: data.table
  
  ### Setup and tests
      
      message("Creating sumtable")
      message(' -- running some initial tests')
      
      dat <- as.data.table(dat)
      exp.cols <- use.cols
      
      if(!is.null(counts)){
        if(is.data.frame(counts)){
          counts <- as.data.table(counts)
        } else {
          stop("counts are not a correctly formatted data.frame or data.table. Either correct the counts, or use counts = NULL")
        }
      }
      
      if(!is.null(parent.col)){
        
        test <- list()
        test2 <- c()
        
        for(i in unique(dat[[pop.col]])){
          # i <- unique(dat[[pop.col]])[1]
          test[[i]] <- unique(unlist(dat[dat[[pop.col]] == i,parent.col, with = FALSE]))
          if(length(test[[i]]) > 1){
            test2[i] <- TRUE
          } else {
            test2[i] <- FALSE
          }
        }
        
        if(any(test2)){
          message('Error: some subpopulations have multiple parent populations listed')
          message(' ')
          for(i in names(test)){
            message(paste0(' -- Parents for ', i, ': '))
            message(paste0('       ', test[[i]]))
          }
          message(' ')
          stop('Either adjust the dataset, or use parent.col = NULL')
        }
        rm(test)
        rm(test2)
      }
  
  ### Test data
      
      # dat <- Spectre::demo.clustered
      # sample.col <- 'Sample'
      # pop.col <- 'FlowSOM_metacluster'
      # annot.cols <- c('Group', 'Batch')
      # exp.cols <- c('CD4_asinh', 'Ly6C_asinh')
      # parent.col <- 'Population'
      # func = 'median'
      # sep <- ' -- '
      # counts <- data.table('A' = unique(dat$Sample),'B' = c(rep(100000, 6), rep(1000000, 6)))
      # perc.pos <- data.table('A' = c('CD4_asinh', 'Ly6C_asinh'), 'B' = c(2.5,2.5))
      # double.pos <- c('CD4_asinh', 'Ly6C_asinh')

  ### Cell proportions
      
      message(" -- calculting cell proportions")
      
      # props <- data.frame(with(dat, table(Population, Sample)))
      props <- data.frame(with(dat, table(dat[[pop.col]], dat[[sample.col]])))
      names(props) <- c(pop.col, sample.col, 'nrows')
      props$NAME <- paste0(props[[pop.col]], ' -- ', props[[sample.col]])
      props <- props[,c(4,1,2,3)]
      
      for(i in unique(props[[sample.col]])){
        props[props[[sample.col]] == i,'Percent of sample'] <- props[props[[sample.col]] == i,'nrows'] / sum(props[props[[sample.col]] == i,'nrows']) * 100
      }
      
      template <- data.table('NAME' = props$NAME)
      template
  
  ### Cell counts
      
      if(!is.null(counts)){
        message(" -- calculating cell counts")
        for(i in unique(props[[sample.col]])){
          props[props[[sample.col]] == i,'Cells per sample'] <- props[props[[sample.col]] == i,'Percent of sample'] * counts[counts[[1]] == i,2][[1]] / 100
        }
      }
      
  ### Proportion adjust
  
      props[[sample.col]] <- NULL
      props[[pop.col]] <- NULL
      
      template <- do.add.cols(template, 'NAME', props, 'NAME', show.status = FALSE)
    
          # template
          # 
          # test <- tidyr::separate(template, 'NAME', sep = ' -- ', into = c(pop.col, sample.col))
          # test$nrows <- NULL
          # test <- reshape(test, idvar = sample.col, timevar = pop.col, direction = "wide", sep = sep)
          # test <- as.data.table(test)
          # test    
      
      gc()
  
  ### Percent of parent
      
      if(!is.null(parent.col)){
        
        message(" -- calculting percent of parent")
        
        res.list <- list()
        
        for(i in unique(dat[[parent.col]])){
          # i <- unique(dat[[parent.col]])[1]
          
          message(paste0("   ... ", i))
          
          rws <- dat[[parent.col]] == i
          
          tmp <- dat[rws,]
          unique(tmp[[pop.col]])
          
          res.1 <- data.frame(with(tmp, table(tmp[[parent.col]], tmp[[sample.col]])))
          res.2 <- data.frame(with(tmp, table(tmp[[pop.col]], tmp[[sample.col]])))
          
          res.2 <- do.filter(res.2, 'Var1', unique(tmp[[pop.col]]))
          
          names(res.1)[3] <- 'ParentCounts'
          names(res.2)[3] <- 'PopCounts'
          
          res.3 <- do.add.cols(res.2, 'Var2', res.1[,2:3], 'Var2', show.status = FALSE)
          res.3$NAME <- paste0(res.3[['Var1']], ' -- ', res.3[['Var2']]) 
          res.3[[paste0('Perc of parent')]] <- res.3$PopCounts / res.3$ParentCounts * 100
          res.3$Var1 <- NULL
          res.3$Var2 <- NULL
          res.3$PopCounts <- NULL
          res.3$ParentCounts <- NULL
          res.list[[i]] <- res.3
          
          rm(res.1)
          rm(res.2)
          rm(res.3)
          rm(tmp)
          rm(rws)
        }
        
        template <- do.add.cols(template, 'NAME', rbindlist(res.list), 'NAME', show.status = FALSE)
        
        rm(res.list)
        gc()
      }
  
  ### Expression
      
      if(!is.null(exp.cols)){
        
        message(" -- calculting expression levels")
        
        if(func == 'median'){
          exps <- dat[, lapply(.SD, median), by = c(sample.col, pop.col), .SDcols = exp.cols]
        }
        if(func == 'mean'){
          exps <- dat[, lapply(.SD, mean), by = c(sample.col, pop.col), .SDcols = exp.cols]
        }
        if(func == 'sum'){
          exps <- dat[, lapply(.SD, sum), by = c(sample.col, pop.col), .SDcols = exp.cols]
        }
        
        names(exps)[-c(1:2)] <- paste0('Exp ', names(exps)[-c(1:2)])
        exps$NAME <- paste0(exps[[pop.col]], ' -- ', exps[[sample.col]])
        exps[[sample.col]] <- NULL
        exps[[pop.col]] <- NULL
        
        template <- do.add.cols(template, 'NAME', exps, 'NAME', show.status = FALSE)
        rm(exps)
        gc()
      }
  
  ### Percent positive
      
      if(!is.null(perc.pos)){
        
        message(" -- calculting percent positive")
        
        alt <- dat[,c(sample.col, pop.col, perc.pos[[1]]), with = FALSE]
        
        for(i in c(1:nrow(perc.pos))){
          # i <- 1
          
          message(paste0("     ... ", perc.pos[[1]][i]))
          
          mrk <- perc.pos[i,1][[1]]
          val <- perc.pos[i,2][[1]]
          
          alt[alt[[mrk]] > val,paste0('POS ', mrk)] <- TRUE
          alt[alt[[mrk]] <= val,paste0('POS ', mrk)] <- FALSE
        }
        
        alt.res <- alt[, lapply(.SD, sum), by = c(sample.col, pop.col), .SDcols = paste0('POS ', perc.pos[[1]])]
        alt.res$NAME <- paste0(alt.res[[pop.col]], ' -- ', alt.res[[sample.col]])
        alt.res[[sample.col]] <- NULL
        alt.res[[pop.col]] <- NULL
        alt.res <- do.add.cols(alt.res, 'NAME', props[,c('NAME', 'nrows')], 'NAME', show.status = FALSE)
        alt.res[,paste0('PROP POS ', perc.pos[[1]])] <- alt.res[,paste0('POS ', perc.pos[[1]]), with = FALSE] / alt.res$nrows * 100
        alt.res <- alt.res[,c('NAME', names(alt.res)[grepl('PROP POS ', names(alt.res))]), with = FALSE]
        
        template <- do.add.cols(template, 'NAME', alt.res, 'NAME', show.status = FALSE)
        rm(alt.res)
        rm(alt)
        gc()
      }
      
  ### Multiple positives
      
      if(!is.null(double.pos)){
        
        message(" -- calculting multiple positive")
        
        alt <- dat[,c(sample.col, pop.col, double.pos), with = FALSE]
        
        for(i in c(1:length(double.pos))){
          # i <- 1
          
          mrk <- double.pos[i]
          val <- perc.pos[perc.pos[[1]] == mrk,2][[1]]
          
          alt[alt[[mrk]] > val,paste0('EXP-', mrk)] <- TRUE
          alt[alt[[mrk]] <= val,paste0('EXP-', mrk)] <- FALSE
        }
        
        pos.cols <- names(alt)[grepl('EXP-', names(alt))]
        
        tmp.list <- list()
        
        if(length(pos.cols) == 2){
          tmp.list[[paste0(pos.cols[1], ' POS', ' + ', pos.cols[2], ' POS')]] <- alt[alt[[pos.cols[1]]] == TRUE &
                                                                                       alt[[pos.cols[2]]] == TRUE,]
          
          tmp.list[[paste0(pos.cols[1], ' POS', ' + ', pos.cols[2], ' NEG')]] <- alt[alt[[pos.cols[1]]] == TRUE &
                                                                                       alt[[pos.cols[2]]] == FALSE,]
          
          tmp.list[[paste0(pos.cols[1], ' NEG', ' + ', pos.cols[2], ' POS')]] <- alt[alt[[pos.cols[1]]] == FALSE &
                                                                                       alt[[pos.cols[2]]] == TRUE,]
          
          tmp.list[[paste0(pos.cols[1], ' NEG', ' + ', pos.cols[2], ' NEG')]] <- alt[alt[[pos.cols[1]]] == FALSE &
                                                                                       alt[[pos.cols[2]]] == FALSE,] 
        }
        
        if(length(pos.cols) == 3){
          tmp.list[[paste0(pos.cols[1], ' POS', ' + ', pos.cols[2], ' POS', ' + ', pos.cols[3], ' POS')]] <- alt[alt[[pos.cols[1]]] == TRUE &
                                                                                                                   alt[[pos.cols[2]]] == TRUE &
                                                                                                                   alt[[pos.cols[3]]] == TRUE,]
          
          tmp.list[[paste0(pos.cols[1], ' POS', ' + ', pos.cols[2], ' POS', ' + ', pos.cols[3], ' NEG')]] <- alt[alt[[pos.cols[1]]] == TRUE &
                                                                                                                   alt[[pos.cols[2]]] == TRUE &
                                                                                                                   alt[[pos.cols[3]]] == FALSE,]
          
          tmp.list[[paste0(pos.cols[1], ' POS', ' + ', pos.cols[2], ' NEG', ' + ', pos.cols[3], ' POS')]] <- alt[alt[[pos.cols[1]]] == TRUE &
                                                                                                                   alt[[pos.cols[2]]] == FALSE &
                                                                                                                   alt[[pos.cols[3]]] == TRUE,]
          
          tmp.list[[paste0(pos.cols[1], ' NEG', ' + ', pos.cols[2], ' POS', ' + ', pos.cols[3], ' POS')]] <- alt[alt[[pos.cols[1]]] == FALSE &
                                                                                                                   alt[[pos.cols[2]]] == TRUE &
                                                                                                                   alt[[pos.cols[3]]] == TRUE,]
          
          tmp.list[[paste0(pos.cols[1], ' POS', ' + ', pos.cols[2], ' NEG', ' + ', pos.cols[3], ' NEG')]] <- alt[alt[[pos.cols[1]]] == TRUE &
                                                                                                                   alt[[pos.cols[2]]] == FALSE &
                                                                                                                   alt[[pos.cols[3]]] == FALSE,]
          
          tmp.list[[paste0(pos.cols[1], ' NEG', ' + ', pos.cols[2], ' POS', ' + ', pos.cols[3], ' NEG')]] <- alt[alt[[pos.cols[1]]] == FALSE &
                                                                                                                   alt[[pos.cols[2]]] == TRUE &
                                                                                                                   alt[[pos.cols[3]]] == FALSE,]
          
          tmp.list[[paste0(pos.cols[1], ' NEG', ' + ', pos.cols[2], ' NEG', ' + ', pos.cols[3], ' POS')]] <- alt[alt[[pos.cols[1]]] == FALSE &
                                                                                                                   alt[[pos.cols[2]]] == FALSE &
                                                                                                                   alt[[pos.cols[3]]] == TRUE,]
          
          tmp.list[[paste0(pos.cols[1], ' NEG', ' + ', pos.cols[2], ' NEG', ' ', pos.cols[3], ' NEG')]] <- alt[alt[[pos.cols[1]]] == FALSE &
                                                                                                                 alt[[pos.cols[2]]] == FALSE &
                                                                                                                 alt[[pos.cols[3]]] == FALSE,]
        }
        
        rm(alt)
        gc()
        
        names(tmp.list)
        
        multi.list <- list()
        
        for(a in names(tmp.list)){
          # a <- names(tmp.list)[1]
          
          message(paste0('    ... ', gsub('EXP-', '', a)))
          
          tp <- tmp.list[[a]]
          tp$MULTIKEY <- TRUE
          multi.res <- tp[, lapply(.SD, sum), by = c(sample.col, pop.col), .SDcols = 'MULTIKEY']
          multi.res$NAME <- paste0(multi.res[[pop.col]], ' -- ', multi.res[[sample.col]])
          multi.res[[sample.col]] <- NULL
          multi.res[[pop.col]] <- NULL
          multi.res <- do.add.cols(multi.res, 'NAME', props[,c('NAME', 'nrows')], 'NAME', show.status = FALSE)
          multi.res[,paste0('PROP MULTIPOS ', a)] <- multi.res[,'MULTIKEY', with = FALSE] / multi.res$nrows * 100
          multi.res$MULTIKEY <- NULL
          multi.res$nrows <- NULL
          
          names(multi.res) <- gsub('EXP-', '', names(multi.res))
          
          template <- do.add.cols(template, 'NAME', multi.res, 'NAME', show.status = FALSE)
          
          rm(tp)
          rm(multi.res)
          gc()
        }
        
        rm(multi.list)
        gc()
      }
  
  ### Wrap up
      
      message(' -- wrapping up')
      
      final <- tidyr::separate(template, 'NAME', sep = ' -- ', into = c(pop.col, sample.col))
      final$nrows <- NULL
      final <- reshape(final, idvar = sample.col, timevar = pop.col, direction = "wide", sep = sep)
      final <- as.data.table(final)
      
      if(!is.null(annot.cols)){
        res.cols <- names(final)[-1]
        
        ann <- data.table()
        
        for(i in unique(dat[[sample.col]])){
          # i <- unique(dat[[sample.col]])[1]
          tp <- dat[dat[[sample.col]] == i,]
          tp <- tp[1,c(sample.col, annot.cols), with = FALSE]  
          ann <- rbind(ann, tp)
        }
        
        final <- do.add.cols(final, sample.col, ann, sample.col, show.status = FALSE)
        final <- final[,c(sample.col, annot.cols, sort(res.cols)), with = FALSE] 
      }
  
  ### Return
  
      message(' -- sumtable complete!')
      return(final)
  
}

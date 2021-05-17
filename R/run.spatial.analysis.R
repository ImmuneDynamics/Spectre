#' run.spatial.analysis
#'
#' @param dat NO DEFAULT. Data.table
#' @param sample.col NO DEFAULT. Column that denotes 'samples'
#' @param pop.col NO DEFAULT. Column that denotes 'samples'
#' @param annot.cols DEFAULT = NULL. Annotation columns.
#' @param region.col DEFAULT = NULL. Create a 'total' by default, and add specific ones if requested here
#' @param area.table DEFAULT = NULL. Calculate 'total' automatically
#' @param adj.dist DEFAULT = 100
#' @param x.col DEFAULT = 'x'
#' @param y.col DEFAULT = 'y'
#' @param distribution DEFAULT = TRUE
#' @param composition DEFAULT = TRUE
#' @param counts DEFAULT = TRUE
#' @param counts.per.area DEFAULT = TRUE
#' @param distance DEFAULT = TRUE
#' @param adjacency DEFAULT = TRUE
#' @param func DEFAULT = 'mean'
#'
#' @import data.table
#'
#' @export

# run.spatial.analysis

run.spatial.analysis <- function(dat,
                                 sample.col,
                                 pop.col,
                                 
                                 annot.cols = NULL,
                                 region.col = NULL, # Create a 'total' by default, and add specific ones if requested here
                                 area.table = NULL, # calculate 'total' automatically
                                 
                                 adj.dist = 100,
                                 
                                 x.col = 'x',
                                 y.col = 'y',
                                 
                                 distribution = TRUE,
                                 composition = TRUE,
                                 counts = TRUE,
                                 counts.per.area = TRUE,
                                 distance = TRUE,
                                 adjacency = TRUE,
  
                                 func = 'mean'){

  ### Test data

      # dat <- cell.dat
      # sample.col <- "ROI"
      # pop.col <- "CellType"
      # region.col <- "Region"
      # 
      # area.table <- area.table
      # adj.dist = 100
      # 
      # annot.cols = 'Group'
      # 
      # distribution = TRUE
      # composition = TRUE
      # counts = TRUE
      # counts.per.area = TRUE
      # distance = TRUE
      # adjacency = TRUE
      # 
      # x.col = 'x'
      # y.col = 'y'
      # 
      # func = 'mean'

  ### Function for NAs

      # ###################### NA to 0 ###################
      do.rmv.na = function(dat) {
        # either of the following for loops

        # by name :
        for (j in names(dat))
          set(dat,which(is.na(dat[[j]])),j,0)

        # or by number (slightly faster than by name) :
        # for (j in seq_len(ncol(dat)))
        #   set(dat,which(is.na(dat[[j]])),j,0)
      }
      # ################################################

  ### Setup

      pops <- unique(dat[[pop.col]])
      samples <- unique(dat[[sample.col]])
      regions <- unique(dat[[region.col]])
      
      setorderv(dat, sample.col)
      setorderv(dat, pop.col)
      
      combs <- gtools::permutations(n = length(pops), r = 2, v = pops, repeats.allowed = TRUE)
      combs # x, y
      
  ### Loop PER SAMPLE

      all.counts <- list()

      for(i in samples){
        # i <- samples[[1]]

        message('Processing ', i)
        
        ### Setup
        
            ## Subset sample data
            samp.dat <- dat[dat[[sample.col]] == i,]
    
            ## Sample front matter
            samp.front <- samp.dat[1,c(sample.col, annot.cols), with = FALSE]
            
            ## Preparation
            samp.counts <- as.data.table(regions)
            names(samp.counts)[1] <- "REGION"
    
            ## Loop for each population across regions
    
            reg.res <- as.data.table(pops)
            names(reg.res)[1] <- "POPULATIONS"
    
            for(a in regions){
              # a <- regions[[1]]
    
              reg.dat <- samp.dat[samp.dat[[region.col]] == a,]
              counts <- reg.dat[, .(count = .N), by = pop.col]
              names(counts)[1] <- "POPULATIONS"
              names(counts)[2] <- a
    
              reg.res <- do.add.cols(reg.res, "POPULATIONS", counts, "POPULATIONS", show.status = FALSE)
              do.rmv.na(reg.res)
            }

        ### DISTRIBUTION -- Type A -- for each cell type, where is it located (row proportions) 
    
            if(distribution == TRUE){
              
              message(' -- Calculating distribution')
              
              a.res <- data.table()
              
              for(o in c(1:nrow(reg.res))){
                # o <- 1
                nme <- reg.res[o,1]
                
                rw.ttl <- sum(reg.res[o,-1])
                res <- reg.res[o,-1]/rw.ttl
                res <- res*100
                
                a.res <- rbind(a.res, cbind(nme, res))
                
                rm(o)
                rm(nme)
                rm(rw.ttl)
                rm(res)
              }
              
              do.rmv.na(a.res)
              a.res.long <- melt(setDT(a.res), id.vars = c("POPULATIONS"), variable.name = "REGION")
              a.res.long
              
              a.res.long.new <- data.table()
              a.res.long.new$measure <- paste0("Distribution of ", a.res.long$POPULATIONS, " in ", a.res.long$REGION, " -- Percent of cell type in sample")
              a.res.long.new$counts <- a.res.long$value
              
              a.res.long.new <- dcast(melt(a.res.long.new, id.vars = "measure"), variable ~ measure)
              
              a.res.long.new$variable <- NULL
              
              # a.res.long.new$variable <- i
              # names(a.res.long.new)[1] <- sample.col
              
            }
            

        ### COMPOSITION --  Type B -- for each region, what cells are in it (column proportions)
    
            if(composition == TRUE){
              
              message(' -- Calculating composition')
              
              b.res <- as.data.table(reg.res[,1])
              
              for(o in c(2:length(names(reg.res)))){
                # o <- 2
                nme <- names(reg.res)[o]
                
                col.ttl <- sum(reg.res[,..o])
                res <- reg.res[,..o]/col.ttl
                res <- res*100
                
                b.res <- cbind(b.res, res)
                
                rm(o)
                rm(nme)
                rm(col.ttl)
                rm(res)
              }
              
              do.rmv.na(b.res)
              b.res.long <- melt(setDT(b.res), id.vars = c("POPULATIONS"), variable.name = "REGION")
              b.res.long
              
              b.res.long.new <- data.table()
              b.res.long.new$measure <- paste0("Composition of ", b.res.long$REGION, " - ", b.res.long$POPULATIONS, " -- Percent of cells in region")
              b.res.long.new$counts <- b.res.long$value
              b.res.long.new
              
              b.res.long.new <- dcast(melt(b.res.long.new, id.vars = "measure"), variable ~ measure)
              
              b.res.long.new$variable <- NULL
              
              #b.res.long.new$variable <- i
              #names(b.res.long.new)[1] <- sample.col
            }

        ### COUNTS
    
            message(' -- Calculating cell counts')
            
            reg.res.long <- melt(setDT(reg.res), id.vars = c("POPULATIONS"), variable.name = "REGION")
            reg.res.long
    
            reg.res.long.new <- data.table()
            reg.res.long.new$measure <- paste0("Cell counts in ", reg.res.long$REGION, " - ", reg.res.long$POPULATIONS, " -- Cells per region")
            reg.res.long.new$counts <- reg.res.long$value
            
            reg.res.long.new <- dcast(melt(reg.res.long.new, id.vars = "measure"), variable ~ measure)
            reg.res.long.new$variable <- NULL
            
            # reg.res.long.new$variable <- i
            # names(reg.res.long.new)[1] <- sample.col

        ### COUNTS / AREA
    
            message(' -- Calculating counts/area')
            
            reg.res.by.area <- reg.res
    
            for(u in regions){
              # u <- regions[[1]]
    
              ar <- area.table[area.table[[sample.col]] == i, u, with = FALSE]
              ar <- ar[[1]]
              ar <- as.numeric(ar)
    
              reg.res.by.area[[u]] <- reg.res.by.area[[u]] / ar * 10000 # per 100 um^2
            }
    
            reg.res.area.long <- melt(setDT(reg.res.by.area), id.vars = c("POPULATIONS"), variable.name = "REGION")
            reg.res.area.long
    
            reg.res.area.long.new <- data.table()
            reg.res.area.long.new$measure <- paste0("Cells per area in ", reg.res.area.long$REGION, " - ", reg.res.area.long$POPULATIONS, " -- Cells per 100 um^2 of region")
            reg.res.area.long.new$counts <- reg.res.area.long$value
            
            reg.res.area.long.new <- dcast(melt(reg.res.area.long.new, id.vars = "measure"), variable ~ measure)
            reg.res.area.long.new$variable <- NULL

        ### DISTANCE and ADJACENCY
            
            message(' -- Calculating distance and adjacency')
            
            comb.dist.list <- list()
            comb.adj.list <- list()
            
            for(a in c(1:nrow(combs))){
              # a <- 1
              
              cell.a <- combs[a,1]
              cell.b <- combs[a,2]
              
              nme <- paste0(cell.a, " - ", cell.b)
              message("   ---- ", paste0(cell.a, " - ", cell.b))
              
              ## Setup
              
              pop.temp <- do.filter(samp.dat, pop.col, c(cell.a, cell.b))
              cell.types <- pop.temp[[pop.col]]
              
              pop.temp <- pop.temp[,c(x.col, y.col),with = FALSE]
              
              A <- pop.temp[which(cell.types == cell.a),]
              B <- pop.temp[which(cell.types == cell.b),]
              
              ## Calculate distance
              
              # https://stackoverflow.com/questions/48117286/distance-between-two-sets-of-points
              res <- proxy::dist(as.data.frame(A), as.data.frame(B))
              
              ## Calculate adjaceny
              
              adj.res <- res < adj.dist
              table(adj.res)
              
              # table(adj.res)["TRUE"] / c(nrow(res) * ncol(res))     # SAME RESULT AS mean(adj.res, na.rm = TRUE)
              
              ## Save results
              
              comb.dist.list[[nme]] <- mean(res)
              comb.adj.list[[nme]] <- mean(adj.res, na.rm = TRUE)
            }

            comb.dist <- as.data.table(comb.dist.list)
            comb.adj <- as.data.table(comb.adj.list)
            
            names(comb.dist) <- paste0('Av. Distance. - ', names(comb.dist))
            names(comb.adj) <- paste0('Av. Num. Neighbours < ', adj.dist, ' - ', names(comb.adj))
            
        ### WRAP UP for sample

            sample.all <- cbind(samp.front, 
                                reg.res.long.new, 
                                reg.res.area.long.new, 
                                a.res.long.new, 
                                b.res.long.new,
                                comb.dist,
                                comb.adj)
            
            all.counts[[i]] <- sample.all
            
            # all.counts <- rbind(all.counts, cbind(reg.res.long.new, 
            #                                       reg.res.area.long.new, 
            #                                       a.res.long.new, 
            #                                       b.res.long.new))
    
            rm(reg.res.long.new)
            rm(reg.res.area.long.new)
            rm(a.res.long.new)
            rm(b.res.long.new)
            rm(comb.dist)
            rm(comb.adj)
            rm(sample.all)
      }
      
  ### Wrap up
      
      final <- rbindlist(all.counts, fill = TRUE)
      return(final)

}

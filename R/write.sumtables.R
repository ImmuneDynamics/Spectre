#' write.sumtables
#'
#' @usage make.sumtables(x, ...)
#'
#' @param x NO DEFAULT. A data.table containing cells (rows) vs features/markers (columns). One column must represent sample names, and another must represent populations/clusters.
#' @param sample.col NO DEFAULT. Character. Name of the sample column (e.g. "Sample").
#' @param pop.col NO DEFAULT. Character. Name of the population/cluster column (e.g. "Population", "Cluster").
#' @param measure.col NO DEFAULT. A character or numeric vector indicating the columns to be measured (e.g. cellular columns -- c("CD45", "CD3e") etc).
#'
#' @param annot.col DEFAULTS TO NULL. A character or numeric vector indicating the columns to be included as annotation columns (e.g. cellular columns -- c("Batch", "Group") etc). If groups are present, this column must be present here also.
#' @param group.col DEFAULTS TO NULL. Character. Which column represent groups (e.g. "Group"). This is for dividing data within the function. If you wish for groups to be included in the annotations, you will need to also include it in the annot.col argument.

#' @param do.frequencies DEFAULTS TO TRUE. Do you wish to create cell proportion and/or cell count results?
#' @param cell.counts DEFAULTS TO NULL. If you wish to generate cell.count results, a vector of cell counts (e.g. c(1000, 1500, 2439,)) representing the cell counts in each of the samples. Must be entered in the order the that unique sample names appear in the dataset.
#' @param do.mfi.per.sample DEFAULTS TO TRUE. Do you wish to generate MFI data (markers vs clusters) for each sample?
#' @param do.mfi.per.marker Currently inactive, set to FALSE.
#'
#' @param perc.pos.markers DEFAULTS to NULL. A vector of column names of calculating percent positive summary stats.
#' @param perc.pos.cutoff DEFAULTS to NULL. A vector of 'positive' cut-off values for the markers defined in perc.pos.markers. Must be in same order.
#'
#' @param mfi.type DEFAULTS TO "median". Can be "median" or "mean". Defines the type of function for calculating MFI data.
#' @param path DEFAULTS TO getwd(). Defines the directory for write CSV files.
#'
#' @return ...
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @examples make.sumtables()
#'
#' @export

write.sumtables <- function(x,
                           sample.col, # name
                           pop.col, # name
                           measure.col, # numbers

                           ## Optional extras -- with implications
                           annot.col = NULL, # numbers # takes the first value from each 'sample' -- as these are annotations that are entered per 'cell', but should be equivalent for each sample
                           group.col = NULL, # name

                           ## Functions
                           do.frequencies = TRUE,
                           cell.counts = NULL,
                           do.mfi.per.sample = TRUE,
                           do.mfi.per.marker = FALSE, # coming soon

                           ## Functions for percentage positive calculations
                           perc.pos.markers = NULL,
                           perc.pos.cutoff = NULL,

                           ## Other defaults
                           mfi.type = "median",
                           path = getwd())
{

#####################################################################
#### Required packages
#####################################################################

  ### Check that necessary packages are installed
      if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
      if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')

  ### Require packages
      require(Spectre)
      require(data.table)

#####################################################################
#### Data for testing
#####################################################################

      # library(Spectre)
      # library(data.table)
      # library(tidyr)
      # x = as.data.table(Spectre::demo.clustered)
      # x <- as.data.frame(x)
      # sample.col = "Sample"
      # pop.col = "FlowSOM_metacluster"
      #
      # measure.col = names(x[c(2,5:6,8:9,11:13,16:19,21:30,32)])
      # annot.col = names(x[c(33:34,36:37)])
      # group.col = "Group"
      # cell.counts = c(rep(2.0e+07, 6), rep(1.8e+07, 6))
      #
      # do.frequencies = TRUE
      # do.mfi.per.sample = TRUE
      # do.mfi.per.marker = TRUE
      #
      # perc.pos.markers = c("BV711.SCA.1","APC.BrdU")
      # perc.pos.cutoff = c(580, 450)
      #
      # mfi.type = "median"
      # path = setwd("/Users/thomasa/Desktop/")

#####################################################################
#### 1. Initial setup
#####################################################################

      ## Order the data by pop and then sample
          # x <- x.start

          message("SumTables - initial setup started")

          x <- as.data.frame(x)

          if(!is.null(cell.counts)){
            count <- data.frame(unique(x[sample.col]), cell.counts)
            count <- count[order(count[1]),]
          }

          x <- x[order(x[pop.col]),]
          x <- x[order(x[sample.col]),]

      ## Unique lists
          all.samples <- unique(x[[sample.col]])
          all.pops <- unique(x[[pop.col]])

          if(!is.null(group.col)){
            all.groups <- unique(x[[group.col]])
          }

          #sample.col.num <- match(sample.col,names(x))
          #group.col.num <- match(group.col,names(x))

      ## Create table of annot columns for append
          x[,annot.col]

          if(length(annot.col) > 1){
            annots <- list()
            for(i in all.samples){
              temp <- x[x[[sample.col]] == i,]
              annots[[i]] <- temp[1,annot.col]
              # add a warning here for if the values aren't identical
            }
            annots <- rbindlist(annots)
            annots
          }

          if(length(annot.col) == 1){
            annots <- list()
            for(i in all.samples){
              temp <- x[x[[sample.col]] == i,]
              annots[[i]] <- temp[1,annot.col]
            }
            annots <- as.data.frame(unlist(annots))
            names(annots) <- annot.col
            annots <- as.data.table(annots)
          }

      ## Some checks

          x[,sample.col]
          x[,group.col]
          x[,measure.col]
          x[,annot.col]

      ## RENAMING
          colnames(x)[which(names(x) == sample.col)] <- "SAMPLE"
          colnames(x)[which(names(x) == pop.col)] <- "POPULATION"

          x[,"SAMPLE"]
          x[,"POPULATION"]

      message("SumTables - initial setup complete")

#####################################################################
#### 2. PROPORTIONS
#####################################################################

      if(do.frequencies == TRUE){

        message("SumTables - frequencies started")

        ## Set up data
        x.freq <- cbind(POPULATION = x[,"POPULATION"], SAMPLE = x[,"SAMPLE"], x[,measure.col])
        x.freq

        all.samples
        all.pops

        ## Count nrows per sample
        samp.list = list()

        for(i in all.samples) {
          temp <- x.freq[x.freq$SAMPLE == i,]
          samp.list[[i]] <- temp
        }

        samp.list
        #samp.list[[7]]
        #nrow(samp.list[[7]])

        new.sample.list = list()

        for(a in all.samples){
          samp.res = data.frame(Cluster = numeric(0), NumCells= numeric(0))

          for(i in all.pops){
            num.cells <- sum(samp.list[[a]]["POPULATION"] == i) # number of cells per cluster
            res <- data.frame(Cluster = i, NumCells = num.cells)
            samp.res <- rbind(samp.res, res)
          }

          new.sample.list[[a]] <- samp.res
          new.sample.list[[a]]["SAMPLE"] <- a
        }

        ## Combine data
        merged.df <- rbindlist(new.sample.list)
        merged.df <- spread(merged.df, Cluster, NumCells)
        merged.df <- as.data.frame(merged.df)

        ## Calculate row totals -- to get total cells per sample

        row.totals <- rowSums(merged.df[c(2:length(names(merged.df)))])
        as.matrix(row.totals)

        per <- apply(merged.df[c(2:length(names(merged.df)))], 1, function(i) i/sum(i)) # https://stackoverflow.com/questions/35691205/divide-each-each-cell-of-large-matrix-by-sum-of-its-row
        per <- t(per)

        per
        per.100 <- per*100

        as.matrix(rowSums(per.100))

        x.freq.res <- cbind("Measurement" = rep("Percentage of total cells", nrow(per.100)), "SAMPLE" = merged.df$SAMPLE, annots, "Row total" = rowSums(per.100), per.100)
        names(x.freq.res)[names(x.freq.res) == "SAMPLE"] <- sample.col
        x.freq.res

        ## Save data
        setwd(path)
        write.csv(x.freq.res, "SumTable-Proportions.csv", row.names=FALSE)

        message("SumTables - frequencies complete")
      }

#####################################################################
#### 3. CELL COUNTS
#####################################################################

  if(do.frequencies == TRUE){
    if(!is.null(cell.counts)){

      message("SumTables - cell counts started")

      ## Checks
      check <- length(cell.counts) == length(unique(all.samples))
      if(check == FALSE){
        message("Please note: you do not have the same number of `cell counts` as you do `samples``, so cell number calculations not performed. Please check your entries and try again.")
      }

      ## Calculations
      cellsper <- sweep(per,MARGIN=1,cell.counts,`*`)

      x.count.res <- cbind("Measurement" = rep("Cells per sample", nrow(cellsper)), "SAMPLE" = merged.df$SAMPLE, annots, "Row total" = rowSums(cellsper), cellsper)
      names(x.count.res)[names(x.count.res) == "SAMPLE"] <- sample.col
      x.count.res

      ## Save data
      setwd(path)
      write.csv(x.count.res, "SumTable-CellCounts.csv", row.names=FALSE)

      message("SumTables - cell counts complete")
    }
  }

#####################################################################
#### 4. MFI PER SAMPLE (each table = cluster x marker)
#####################################################################

    if(do.mfi.per.sample == TRUE){

      message("SumTables - MFI per sample started")
      ## All cells

      x.mfi <- cbind(POPULATION = x[,"POPULATION"], x[,measure.col])
      x.mfi

      data.MFI <- aggregate(. ~POPULATION, data = x.mfi, FUN=mfi.type)
      colnames(data.MFI)[1] <- pop.col

      setwd(path)
      dir.create("SumTable-MFI-PerSample", showWarnings = FALSE)
      setwd("SumTable-MFI-PerSample")
      write.csv(data.MFI, "SumTable-MFI-AllSamples.csv", row.names=FALSE)

      ## Per sample

      for(i in all.samples){
        #i <- "01_Mock_01"
        sample.temp <- x[x[["SAMPLE"]] == i,]

        x.mfi.per.sample <- cbind(POPULATION = sample.temp[,"POPULATION"], sample.temp[,measure.col])
        x.mfi.per.sample

        data.MFI.per.sample <- aggregate(. ~POPULATION, data = x.mfi.per.sample, FUN=mfi.type)
        colnames(data.MFI.per.sample)[1] <- pop.col

        setwd(path)
        dir.create("SumTable-MFI-PerSample", showWarnings = FALSE)
        setwd("SumTable-MFI-PerSample")
        write.csv(data.MFI.per.sample, paste0("SumTable-MFI-Sample-", i, ".csv"), row.names=FALSE)
      }

      ## Per group

      if(!is.null(group.col)){
        for(i in all.groups){
          #i <- "Mock"
          grp.temp <- x[x[[group.col]] == i,]

          x.mfi.per.group <- cbind(POPULATION = grp.temp[,"POPULATION"], grp.temp[,measure.col])
          x.mfi.per.group

          data.MFI.per.group <- aggregate(. ~POPULATION, data = x.mfi.per.group, FUN=mfi.type)
          colnames(data.MFI.per.group)[1] <- pop.col

          setwd(path)
          dir.create("SumTable-MFI-PerSample", showWarnings = FALSE)
          setwd("SumTable-MFI-PerSample")
          write.csv(data.MFI.per.group, paste0("SumTable-MFI-Group-", i, ".csv"), row.names=FALSE)
        }
      }
      setwd(path)
      message("SumTables - MFI per sample complete")
    }

#####################################################################
#### 5. MFI PER MARKER (each table = cluster x sample)
#####################################################################

    if(do.mfi.per.marker == TRUE){

      message("SumTables - MFI per marker started")

    ### Setup

        ## Define marker names
        markers <- names(x[,measure.col])
        markers <- sort(markers, decreasing = FALSE)

        ## Define column names
        populations <- unique(x$POPULATION)
        populations <- sort(populations, decreasing = FALSE)

        ## Define data
        temp <- cbind("SAMPLE" = x[,"SAMPLE"], "POPULATION" = x[,"POPULATION"], x[,measure.col])

    ### Loop for each marker

        for(i in markers){
          mkr <- data.frame(row.names = populations)
          test <- data.table("SAMPLE" = x[,"SAMPLE"],
                             "POPULATION" = x[,"POPULATION"],
                             "MARKER" = x[,i])

          ## Sample loop
          for(a in all.samples){

            # u <- test[,c(2,3)]
            #
            # ad <- aggregate(. ~POPULATION, data = u, FUN=mfi.type)
            #
            # xp <- test[test[["SAMPLE"]] == a,]
            # xp <- xp[,c(2,3)] # don't need to see the 'sample

            xp <- test[test[["SAMPLE"]] == a,]
            xp <- xp[,c(2,3)] # don't need to see the 'sample
            ad <- aggregate(. ~POPULATION, data = xp, FUN=mfi.type)

            ## Modification to adjust for any missing clusters (i.e. a sample where cluster 12 is missing, etc)

            if(length(populations) != nrow(ad)){

              strt <- as.data.frame(x = populations)
              names(strt) <- "POPULATION"

              ad <- merge(strt,ad,all = TRUE)
            }

            names(ad)[2] <- as.character(a)
            mkr <- cbind(mkr, ad[2])
          }

          ## Final prep
          mkr <- t(mkr)
          mkr <- as.data.frame(mkr)

          mkr.res <- cbind("Measurement" = rep(paste0("MFI of ", i), nrow(mkr)), "SAMPLE" = rownames(mkr), annots, "Blank" = rep("Blank", nrow(mkr)),mkr)
          mkr.res

          names(mkr.res)[2] <- sample.col
          mkr.res

          ## Save MARKER CSV
          setwd(path)
          write.csv(mkr.res, paste0("SumTable-MFI-", i, ".csv"), row.names = FALSE)
        }

        message("SumTables - MFI per marker complete")
        setwd(path)
    }

#####################################################################
#### 6. Percent positive
#####################################################################

    if(!is.null(perc.pos.markers)){
      if(!is.null(perc.pos.cutoff)){
        message("SumTables - percent positive started")

        setwd(path)

        ### Setup
            x <- as.data.table(x)

            clusters <- all.pops
            clusters

            samples <- all.samples
            samples

            MarkerCutoffs <- data.frame(perc.pos.markers, perc.pos.cutoff)
            MarkerCutoffs

        ### Loops

          for(i in c(1:nrow(MarkerCutoffs))){ ## EACH CSV IS A MARKER
            # i <- 1
            active.marker <- MarkerCutoffs[i,1] #,1 is name, ,2 is cutoff
            active.marker <- as.character(active.marker)
            active.marker

            active.marker.cutoff <- MarkerCutoffs[i,2]
            active.marker.cutoff

            df <- data.frame(matrix(ncol = (1 + length(clusters)), nrow = 0))
            colnames(df) <- c("SAMPLE", as.vector(clusters))
            df

            for(a in samples){ ## EACH ROW IS A SAMPLE
              #a <- "01_Mock_01"
              active.sample <- a

              x.samp <- x[x[["SAMPLE"]] == a,]
              x.samp

              all.cluster.res <- list()

              for(o in clusters){  ## EACH COLUMN IS THE RESULTS OF A CLUSTER
                #o <- 1
                active.cluster <- o

                active.sample.cluster <- x.samp[x.samp[["POPULATION"]] == o,] # Sample == a, Cluster == o
                total <- nrow(active.sample.cluster) # Number of total events of this cluster from this sample

                pos <- active.sample.cluster[active.sample.cluster[[active.marker]] >active.marker.cutoff ,]
                num.pos <- nrow(pos)

                freq <- num.pos / total
                freq <- freq*100
                freq

                res <- data.frame(Cluster = o, Freq = freq)
                res

                all.cluster.res <- rbind(all.cluster.res, res)
              }

              active.sample
              all.cluster.res

              df[nrow(df) + 1,] = c(active.sample, as.vector(all.cluster.res[,2]))
              df

            }

          ## Save
          nms <- df$SAMPLE

          df$SAMPLE <- NULL
          df <- cbind("Measurement" = rep(paste0("Percent ", active.marker, " positive"), nrow(df)),
                      "Sample" = nms,
                      annots,
                      "Blank" = rep("Blank", nrow(df)),
                      df)

          #names(df)[2] <- sample.col

          setwd(path)
          write.csv(df, paste0("SumTable-PercPos-", active.marker, ".csv"), row.names = FALSE)

          message("SumTables - ", active.marker, " percent positive - CSV written to disk")
          setwd(path)
         }
      }
    }

message("SumTables complete")

}


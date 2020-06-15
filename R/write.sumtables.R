#' write.sumtables - Export generated data as a .csv file
#'
#' This function summarises generated results and exports them as a .csv file.
#' It can include the calculation of median fluorescence intensity (or equivalent), proportion (%) or cell counts for clusters/subsets of interest.
#' Makes use of the packages 'data.table' and 'tidyr' to handle the data.
#'
#' @param dat NO DEFAULT. A data.table containing cells (rows) vs features/markers (columns). One column must represent sample names, and another must represent populations/clusters.
#' @param sample.col NO DEFAULT. Character. Name of the sample column (e.g. "Sample").
#' @param pop.col NO DEFAULT. Character. Name of the population/cluster column (e.g. "Population", "Cluster").
#' @param measure.col NO DEFAULT. A character or numeric vector indicating the columns to be measured (e.g. cellular columns -- c("CD45", "CD3e") etc).
#'
#' @param annot.col DEFAULT = NULL. A character or numeric vector indicating the columns to be included as annotation columns (e.g. cellular columns -- c("Batch", "Group") etc). If groups are present, this column must be present here also.
#' @param group.col DEFAULT = NULL. Character. Which column represent groups (e.g. "Group"). This is for dividing data within the function. If you wish for groups to be included in the annotations, you will need to also include it in the annot.col argument.

#' @param do.proportions DEFAULT = TRUE. Do you wish to create cell proportion and/or cell count results?
#' @param cell.counts DEFAULT = NULL. If you wish to generate cell.count results, a vector of cell counts (e.g. c(1000, 1500, 2439,)) representing the cell counts in each of the samples. Must be entered in the order the that unique sample names appear in the dataset.
#' @param do.mfi.per.sample DEFAULT = TRUE. Do you wish to generate MFI data (markers vs clusters) for each sample?
#' @param do.mfi.per.marker Currently inactive, set to FALSE.
#'
#' @param perc.pos.markers DEFAULT = NULL. A vector of column names of calculating percent positive summary stats.
#' @param perc.pos.cutoff DEFAULT = NULL. A vector of 'positive' cut-off values for the markers defined in perc.pos.markers. Must be in same order.
#'
#' @param mfi.type DEFAULT = "median". Can be "median" or "mean". Defines the type of function for calculating MFI data.
#' @param path DEFAULT = getwd(). Defines the directory for write CSV files.
#' 
#' @usage make.sumtables(dat, sample.col, pop.col, measure.col, annot.col, group.col, do.proportions, cell.counts, do.mfi.per.sample, do.mfi.per.marker, perc.pos.markers, perc.pos.cutoff, mfi.type, path)
#'
#' @examples
#' # Calculate and export results from demonstration data
#' Spectre::write.sumtables(dat = Spectre::demo.clustered,
#'                          sample.col = "Sample",
#'                          pop.col = "FlowSOM_metacluster",
#'                          measure.col = c(2,5:6,8:9,11:13,16:19,21:30,32),
#'                          annot.col = c(33:34,36:37),
#'                          group.col = "Group",
#'                          cell.counts = c(rep(2.0e+07, 6), rep(1.8e+07, 6)),
#'                          do.mfi.per.marker = TRUE,
#'                          perc.pos.markers = c("BV711.SCA.1","APC.BrdU"),
#'                          perc.pos.cutoff = c(580, 450)
#'                          )
#'
#' @author
#' Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#' @export

write.sumtables <- function(dat,
                           sample.col, # name
                           pop.col, # name
                           measure.col, # numbers

                           ## Optional extras -- with implications
                           annot.col = NULL, # numbers # takes the first value from each 'sample' -- as these are annotations that are entered per 'cell', but should be equivalent for each sample
                           group.col = NULL, # name

                           ## Functions
                           do.proportions = TRUE,
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
      if(!is.element('tidyr', installed.packages()[,1])) stop('tidyr is required but not installed')
  
  ### Require packages
      require(Spectre)
      require(data.table)
      require(tidyr)

#####################################################################
#### Data for testing
#####################################################################

      # library(Spectre)
      # library(data.table)
      # library(tidyr)
      # dat = as.data.table(Spectre::demo.clustered)
      # dat <- as.data.frame(dat)
      # measure.col = names(dat[c(2,5:6,8:9,11:13,16:19,21:30,32)])
      # annot.col = names(dat[c(33:34,36:37)])
      #
      # do.mfi.per.marker = TRUE

#####################################################################
#### 1. Initial setup
#####################################################################

      ## Order the data by pop and then sample
          # dat <- dat.start

          message("SumTables - initial setup started")

          dat <- as.data.frame(dat)

          if(!is.null(cell.counts)){
            count <- data.frame(unique(dat[sample.col]), cell.counts)
            count <- count[order(count[1]),]
          }

          dat <- dat[order(dat[pop.col]),]
          dat <- dat[order(dat[sample.col]),]

      ## Unique lists
          all.samples <- unique(dat[[sample.col]])
          all.pops <- unique(dat[[pop.col]])

          if(!is.null(group.col)){
            all.groups <- unique(dat[[group.col]])
          }

          #sample.col.num <- match(sample.col,names(dat))
          #group.col.num <- match(group.col,names(dat))

      ## Create table of annot columns for append
          dat[,annot.col]

          if(length(annot.col) > 1){
            annots <- list()
            for(i in all.samples){
              temp <- dat[dat[[sample.col]] == i,]
              annots[[i]] <- temp[1,annot.col]
              # add a warning here for if the values aren't identical
            }
            annots <- data.table::rbindlist(annots)
            annots
          }

          if(length(annot.col) == 1){
            annots <- list()
            for(i in all.samples){
              temp <- dat[dat[[sample.col]] == i,]
              annots[[i]] <- temp[1,annot.col]
            }
            annots <- as.data.frame(unlist(annots))
            names(annots) <- annot.col
            annots <- as.data.table(annots)
          }

      ## Some checks

          dat[,sample.col]
          dat[,group.col]
          dat[,measure.col]
          dat[,annot.col]

      ## RENAMING
          colnames(dat)[which(names(dat) == sample.col)] <- "SAMPLE"
          colnames(dat)[which(names(dat) == pop.col)] <- "POPULATION"

          dat[,"SAMPLE"]
          dat[,"POPULATION"]

      message("SumTables - initial setup complete")

#####################################################################
#### 2. PROPORTIONS
#####################################################################

      if(do.proportions == TRUE){

        message("SumTables - proportions started")

        ## Set up data
        dat.prop <- cbind(POPULATION = dat[,"POPULATION"], SAMPLE = dat[,"SAMPLE"], dat[,measure.col])
        dat.prop

        all.samples
        all.pops

        ## Count nrows per sample
        samp.list = list()

        for(i in all.samples) {
          temp <- dat.prop[dat.prop$SAMPLE == i,]
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
        merged.df <- data.table::rbindlist(new.sample.list)
        merged.df <- tidyr::spread(merged.df, Cluster, NumCells)
        merged.df <- as.data.frame(merged.df)

        ## Calculate row totals -- to get total cells per sample

        row.totals <- rowSums(merged.df[c(2:length(names(merged.df)))])
        as.matrix(row.totals)

        per <- apply(merged.df[c(2:length(names(merged.df)))], 1, function(i) i/sum(i)) # https://stackoverflow.com/questions/35691205/divide-each-each-cell-of-large-matrix-by-sum-of-its-row
        per <- t(per)

        per
        per.100 <- per*100

        as.matrix(rowSums(per.100))

        dat.prop.res <- cbind("Measurement" = rep("Percentage of total cells", nrow(per.100)), "SAMPLE" = merged.df$SAMPLE, annots, "Row total" = rowSums(per.100), per.100)
        names(dat.prop.res)[names(dat.prop.res) == "SAMPLE"] <- sample.col
        dat.prop.res

        ## Save data
        setwd(path)
        data.table::fwrite(dat.prop.res, "SumTable-Proportions.csv", row.names=FALSE)

        message("SumTables - proportions complete")
      }

#####################################################################
#### 3. CELL COUNTS
#####################################################################

  if(do.proportions == TRUE){
    if(!is.null(cell.counts)){

      message("SumTables - cell counts started")

      ## Checks
      check <- length(cell.counts) == length(unique(all.samples))
      if(check == FALSE){
        message("Please note: you do not have the same number of `cell counts` as you do `samples``, so cell number calculations not performed. Please check your entries and try again.")
      }

      ## Calculations
      cellsper <- sweep(per,MARGIN=1,cell.counts,`*`)

      dat.count.res <- cbind("Measurement" = rep("Cells per sample", nrow(cellsper)), "SAMPLE" = merged.df$SAMPLE, annots, "Row total" = rowSums(cellsper), cellsper)
      names(dat.count.res)[names(dat.count.res) == "SAMPLE"] <- sample.col
      dat.count.res

      ## Save data
      setwd(path)
      data.table::fwrite(dat.count.res, "SumTable-CellCounts.csv", row.names=FALSE)

      message("SumTables - cell counts complete")
    }
  }

#####################################################################
#### 4. MFI PER SAMPLE (each table = cluster x marker)
#####################################################################

    if(do.mfi.per.sample == TRUE){

      message("SumTables - MFI per sample started")
      ## All cells

      dat.mfi <- cbind(POPULATION = dat[,"POPULATION"], dat[,measure.col])
      dat.mfi

      data.MFI <- aggregate(. ~POPULATION, data = dat.mfi, FUN=mfi.type)
      colnames(data.MFI)[1] <- pop.col

      setwd(path)
      dir.create("SumTable-MFI-PerSample", showWarnings = FALSE)
      setwd("SumTable-MFI-PerSample")
      data.table::fwrite(data.MFI, "SumTable-MFI-AllSamples.csv", row.names=FALSE)

      ## Per sample

      for(i in all.samples){
        #i <- "01_Mock_01"
        sample.temp <- dat[dat[["SAMPLE"]] == i,]

        dat.mfi.per.sample <- cbind(POPULATION = sample.temp[,"POPULATION"], sample.temp[,measure.col])
        dat.mfi.per.sample

        data.MFI.per.sample <- aggregate(. ~POPULATION, data = dat.mfi.per.sample, FUN=mfi.type)
        colnames(data.MFI.per.sample)[1] <- pop.col

        setwd(path)
        dir.create("SumTable-MFI-PerSample", showWarnings = FALSE)
        setwd("SumTable-MFI-PerSample")
        data.table::fwrite(data.MFI.per.sample, paste0("SumTable-MFI-Sample-", i, ".csv"), row.names=FALSE)
      }

      ## Per group

      if(!is.null(group.col)){
        for(i in all.groups){
          #i <- "Mock"
          grp.temp <- dat[dat[[group.col]] == i,]

          dat.mfi.per.group <- cbind(POPULATION = grp.temp[,"POPULATION"], grp.temp[,measure.col])
          dat.mfi.per.group

          data.MFI.per.group <- aggregate(. ~POPULATION, data = dat.mfi.per.group, FUN=mfi.type)
          colnames(data.MFI.per.group)[1] <- pop.col

          setwd(path)
          dir.create("SumTable-MFI-PerSample", showWarnings = FALSE)
          setwd("SumTable-MFI-PerSample")
          data.table::fwrite(data.MFI.per.group, paste0("SumTable-MFI-Group-", i, ".csv"), row.names=FALSE)
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
        markers <- names(dat[,measure.col])
        markers <- sort(markers, decreasing = FALSE)

        ## Define column names
        populations <- unique(dat$POPULATION)
        populations <- sort(populations, decreasing = FALSE)

        ## Define data
        temp <- cbind("SAMPLE" = dat[,"SAMPLE"], "POPULATION" = dat[,"POPULATION"], dat[,measure.col])

    ### Loop for each marker

        for(i in markers){
          mkr <- data.frame(row.names = populations)
          test <- data.table::data.table("SAMPLE" = dat[,"SAMPLE"],
                             "POPULATION" = dat[,"POPULATION"],
                             "MARKER" = dat[,i])

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
            dat <- data.table::as.data.table(dat)

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

              dat.samp <- dat[dat[["SAMPLE"]] == a,]
              dat.samp

              all.cluster.res <- list()

              for(o in clusters){  ## EACH COLUMN IS THE RESULTS OF A CLUSTER
                #o <- 1
                active.cluster <- o

                active.sample.cluster <- dat.samp[dat.samp[["POPULATION"]] == o,] # Sample == a, Cluster == o
                total <- nrow(active.sample.cluster) # Number of total events of this cluster from this sample

                pos <- active.sample.cluster[active.sample.cluster[[active.marker]] >active.marker.cutoff ,]
                num.pos <- nrow(pos)

                prop <- num.pos / total
                prop <- prop*100
                prop

                res <- data.frame(Cluster = o, Prop = prop)
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


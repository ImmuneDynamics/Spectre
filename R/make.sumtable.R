# make.sumtable

make.sumtable <- function(## Main entries
                           x, # data
                           type, # What kind of sumtables is being generated # Can be: frequencies, expression.per.sample, expression.per.marker

                           ##
                           sample.col, # Samples
                           group.col = NULL, # specify group name (if groups are present)
                           clust.col, # Clusters/populations
                           annot.col.nums = NULL, # Non-cellular markers

                           ## Frequencies
                           cells.per.tissue = NULL, # Calculate cell numbers per cluster/sample/tissue        # Turn OFF if you don't have any cell counts

                           ## Expression
                           fun.type = "median" # Choose summary function for expression levels. Default = "mean". Can be "mean" or "median".

                           ## Groups and fold-change
                           #do.foldchange = FALSE, # Calculate cell number/proportion changes as fold-change         # Turn OFF if you don't have group keywords
                           #ctrl.group = NULL
                           )
{

####################################################################################################################
#### Test data
####################################################################################################################

          # x = cell.dat
          # type = "frequencies"  # Can be: frequencies, expression.per.sample, expression.per.marker
          #
          # sample.col = "Sample"
          # group.name = "Group" # specify group name (if groups are present)
          # clust.col = "FlowSOM_metacluster"
          # annot.col.nums = c(1,33:39) # specify all columns that represent annotations, and NOT CELLULAR PARAMETERS
          #
          # cells.per.tissue = c(rep(2.0e+07, 6), rep(1.8e+07, 6))
          #
          # fun.type = "median"
          #
          # do.foldchange = TRUE # Calculate cell number/proportion changes as fold-change         # Turn OFF if you don't have group keywords
          # ctrl.group = "Mock"

####################################################################################################################
#### Setup
####################################################################################################################

        Starting.dir <- getwd()

        all.sample.names <- as.matrix(unique(x[sample.col]))
        if(is.null(group.col) == FALSE){all.group.names <- as.matrix(unique(x[group.col]))}

        clust.col.num <- which(colnames(x) == clust.col)
        annot.col.nums

        # If cluster column is included in annotation column
        if(any(grepl(clust.col.num, annot.col.nums) == TRUE)){
          annot.col.nums <- setdiff(annot.col.nums, clust.col.num)
          # Now the clustering column has been removed from the annotation columns
        }

        clust.col.num
        annot.col.nums

####################################################################################################################
#### Frequencies
####################################################################################################################

        if(type == "frequencies"){

        ########################## Cells per file ##########################

                  setwd(Starting.dir)
                  dir.create("Output_CellNums", showWarnings = FALSE)

                  data.cellnum <- x  # Specify which column denotes the sample names, and check sample names -- 'table' function will show how many rows (cells) exist per sample
                  data.cellnum <- data.cellnum[order(data.cellnum[sample.col]),] # re-order by sample name (alphabetically)

                  # subset data to cluster number/name, sample number/name, and other parameters that are helpful (e.g. group name, sample name, etc) -- can remove all marker columns (e.g. 141Nd-Ly6G)
                  data.pops <- data.cellnum[c(annot.col.nums, clust.col.num)]
                  samp.name <- data.pops[sample.col]

                  # Split data into separate 'samples'  - specific sample name, all clusters
                  samp.list = list()

                  for(i in all.sample.names) {
                    temp <- subset(data.pops, samp.name == i)
                    samp.list[[i]] <- temp
                  }

                  # some checks
                  samp.list[[1]]
                  nrow(samp.list[[1]])

                  all.samps <- rbindlist(samp.list)
                  all.samps

                  ClustName <- all.samps[[clust.col]]

                  # number of cells per cluster in each sample and order
                  length(as.matrix(unique(ClustName))) # find number of clusters

                  num.clusters <- as.matrix(unique(ClustName)) # list of cluster numbers or names
                  num.clusters <- sort(num.clusters, decreasing = FALSE)

                  new.sample.list = list()

                  for(a in all.sample.names){
                    samp.res = data.frame(Cluster = numeric(0), NumCells= numeric(0), Group=character(0)) # empty dataframe with names

                    for(i in num.clusters){
                      num.cells <- sum(samp.list[[a]][clust.col] == i) # number of cells per cluster
                      #group.label <- samp.list[[a]]$GroupName[1] # picks the group name from the first group of that sample

                      if(is.null(group.col) == FALSE){
                        group.label <- samp.list[[a]][group.col][1,1] # picks the group name from the first group of that sample
                        res <- data.frame(Cluster = i, NumCells = num.cells, Group = group.label) # added group label
                      }

                      if(is.null(group.col) == TRUE){
                        res <- data.frame(Cluster = i, NumCells = num.cells) # added group label
                      }

                      samp.res <- rbind(samp.res, res)
                    }

                    new.sample.list[[a]] <- samp.res
                    new.sample.list[[a]][sample.col] <- a
                  }

                  merged.df <- rbindlist(new.sample.list) ## Concatenate into one large data frame -- what if they have a column conflict??

                  ## Arrange table (long to wide)
                  merged.df <- spread(merged.df, Cluster, NumCells)
                  merged.df <- as.data.frame(merged.df)

                  rownames(merged.df)
                  colnames(merged.df)
                  merged.df[sample.col] # determine which column has the desired 'row names' -- might be first or second column, depending on whether group ID was added

                  ## Save table
                  setwd(Starting.dir)
                  write.csv(merged.df, "Output_CellNums/SumTable_CellsPerFile.csv")

                  # Resulting table represents samples (rows) and clusters (columns). Numbers in each cell represent the number of cells from each cluster in each sample
                  ## Data can now be examined


        ##################################### Cells proportions #####################################

                  prop.df <- merged.df

                  # Save group and
                  if(is.null(group.col) == FALSE){ # Using groups
                    labs <- prop.df[group.col]
                    labs[sample.col] <-  prop.df[sample.col]
                    prop.df[sample.col] <- NULL
                    prop.df[group.col] <- NULL
                  }

                  if(is.null(group.col) == TRUE){ # Not using groups
                    labs <-  prop.df[sample.col]
                    prop.df[sample.col] <- NULL
                  }

                  # create ROW totals (i.e. total cells per sample)
                  row.totals <- rowSums(prop.df)
                  as.matrix(row.totals) # sum(row.totals)

                  cells.per.file <- cbind(labs, row.totals)
                  cells.per.file

                  per <- apply(prop.df, 1, function(i) i/sum(i)) # https://stackoverflow.com/questions/35691205/divide-each-each-cell-of-large-matrix-by-sum-of-its-row
                  per <- t(per)                                                       # 1 = rows, 2 = columns
                  sumtables.percentages <- cbind(labs, per)

                  setwd(Starting.dir)
                  write.csv(sumtables.percentages, "Output_CellNums/SumTable_CellProportions.csv")

        ##################################### Cell counts #####################################


                  if(is.null(cells.per.tissue) == FALSE){

                    check <- length(cells.per.tissue) == length(unique(all.sample.names))
                    if(check == FALSE){
                      print("Error: you do not have the same number of cell counts as samples, cell number calculations not performed")
                    }

                    # times by cell count # input string of 'cell counts'
                    as.matrix(cells.per.tissue)

                    cellsper <- sweep(per,MARGIN=1,cells.per.tissue,`*`)

                    if(is.null(group.col) == FALSE){ # using groups
                      sumtables.cellspertissue <- cbind(labs, cellsper)
                    }
                    if(is.null(group.col) == TRUE){ # NOT using groups
                      sumtables.cellspertissue <- cbind(labs, cellsper)
                    }

                    setwd(Starting.dir)
                    write.csv(sumtables.cellspertissue, "Output_CellNums/SumTable_CellsPerTissue.csv")
                  }


        }


####################################################################################################################
#### MFI cluster x marker (per sample)
####################################################################################################################

      if(type == "expression.per.sample"){

        setwd(Starting.dir)
        dir.create("Output_MFI_per_sample", showWarnings = FALSE)

        ### 5.1 - MFI SUMMARY FOR WHOLE DATASET # Use aggregate to determine MFI of each column, by CLUSTER
            temp <- x

            colnames(temp)[which(names(temp) == clust.col)] <- "CLUSTER"
            fun.type <- noquote(fun.type)

            temp <- temp[-c(annot.col.nums)]
            head(temp)

            data.MFI <- aggregate(. ~CLUSTER, data = temp, FUN=fun.type)
            colnames(data.MFI)[1] <- clust.col

            write.csv(data.MFI, "Output_MFI_per_sample/SumTable_MFI_cluster_x_marker-all_data.csv")


        ## 5.2 - LOOP TO REPEAT FOR EACH SAMPLE
            for(a in all.sample.names){
              temp <- x

              data.subset <- subset(temp, temp[[sample.col]] == a)
              data.subset <- data.subset[-c(annot.col.nums)]
              colnames(data.subset)[which(names(data.subset) == clust.col)] <- "CLUSTER"

              data.subset.MFI <- aggregate(. ~CLUSTER, data = data.subset, FUN=fun.type)
              colnames(data.subset.MFI)[1] <- clust.col # naming the first column

              write.csv(data.subset.MFI, paste("Output_MFI_per_sample/SumTable_MFI_cluster_x_marker-Sample_", a, ".csv", sep=""))
            }

        ## 5.3 - LOOP TO REPEAT FOR EACH GROUP
            if(is.null(group.col) == FALSE){
              for(a in all.group.names){
                temp <- x

                data.subset <- subset(temp, temp[[group.col]] == a)
                data.subset <- data.subset[-c(annot.col.nums)]
                colnames(data.subset)[which(names(data.subset) == clust.col)] <- "CLUSTER"

                data.subset.MFI <- aggregate(. ~CLUSTER, data = data.subset, FUN=fun.type)
                colnames(data.subset.MFI)[1] <- clust.col  # naming the first column

                write.csv(data.subset.MFI, paste("Output_MFI_per_sample/SumTable_MFI_cluster_x_marker-Group_", a, ".csv", sep=""))
              }
            }

      }

####################################################################################################################
#### MFI cluster x sample (per marker)
####################################################################################################################

      if(type == "expression.per.marker"){

        setwd(Starting.dir)
        dir.create("Output_MFI_per_marker", showWarnings = FALSE)

        ### 6.1 - MFI SUMMARY FOR WHOLE DATASET # Use aggregate to determine MFI of each column, by CLUSTER
        temp <- x

        markers <- names(temp)[-c(annot.col.nums, clust.col.num)]
        #colnames(temp)[which(names(temp) == clust.col)] <- "CLUSTER"

        num.clusters <- as.matrix(unique(temp[clust.col])) # find number of clusters
        num.clusters <- sort(num.clusters, decreasing = FALSE)


        ### 6.2 - Cycle through markers

        for (i in c(1:length(unique(markers)))){
          #i <- 1
          z <- markers[i]
          temp[[z]]
          list.per.marker <- list()

          ### 6.3 - Cycle through samples

          for (a in c(all.sample.names)){
            #a <- "01_Mock_01"
            data.subset <- subset(temp, temp[[sample.col]] == a)
            data.subset <- data.subset[-c(annot.col.nums)]

            #data.subset <- aggregate(. ~ CLUSTER, data = data.subset, FUN=fun.type)  # works without this line, but doesn't work with this line if some samples are missing clusters
            ## add something here to deal with missing cluster value. c(1:length_of_clusters) --
            ## if the number of row doesnt not equal the expected number of clusters, add a row for which ever

            ## Replacement for aggregate:
            median.per.cluster <- list()

            for(c in num.clusters){
              #c<-1
              r <- data.subset[data.subset[clust.col] == c,] # selects rows that belong to cluster 'x'
              r
              #r <- r[-c(1:11), ]
              if(fun.type == "mean"){r <- colMeans(r)} # for mean
              if(fun.type == "median"){r <- apply(r, 2, FUN = median)} # for median
              median.per.cluster[[c]] <- r
            }

            median.per.cluster

            #TEST <- rbindlist(list(median.per.cluster), use.names=TRUE)    # doesn't work, don't use
            TEST <- do.call(rbind, unname(median.per.cluster))
            rownames(TEST) <- TEST[,clust.col]
            TEST <- TEST[,-c(which( colnames(TEST)==clust.col))]
            TEST
            data.sub <- as.data.frame(TEST)

            data.sub[[i]]
                #data.subset[i]

            ## Have subseted a marker from one sample
            list.per.marker[[a]] <- data.sub[[i]]  ### This is where the cluster names get replaced by numbers?
                #list.per.marker[[a]] <- data.subset[i]
            list.per.marker
            }

                          #ORIGINAL# lst <- rbindlist(list(list.per.marker)) ## rbind list won't work with samples that are missing some clusters -- tried to use rbind in this case, didn't help

                          #str(list.per.marker)

                          #for(i in c(1:length(list.per.marker))){
                         #   list.per.marker[[i]] <- as.data.frame(list.per.marker[[i]])
                          #}

                          #u <- t(list.per.marker)

                         # rbindlist(list.per.marker)

                          #n <- rbindlist(list(list.per.marker))
                          #n

                          #n <- plyr::rbind.fill(list.per.marker)
                          #n

                          #lst <- rbind(list(list.per.marker))
                          #lst <- as.data.frame(t(lst))
                          #lst

          lst <- as.data.frame(do.call(rbind,list.per.marker))
          lst
          colnames(lst) <- num.clusters
          lst

          write.csv(lst, paste0("Output_MFI_per_marker/SumTable_MFI_per_marker", "_", z, ".csv"))
        }
      }


####################################################################################################################
#### Fold-change
####################################################################################################################


        # ### 4.4 - Calculate FOLD CHANGE (PROPORTIONS)
        #
        # if(do.FOLDCHANGE == 1){
        #   if(use.groups == 1){
        #     sumtables.percentages
        #
        #     ctrl.grp <- subset(sumtables.percentages, Group == ctrls)
        #     ctrl.grp
        #
        #     ctrl.grp.means <- colMeans(ctrl.grp[-c(1:2)])
        #     as.matrix(ctrl.grp.means)
        #
        #     # As above, divide each 'cell' to the vector (by column...)
        #     fold.raw <- t(t(per) / ctrl.grp.means)
        #     fold.raw
        #
        #     sumtables.fold.raw <- cbind(merged.df[c(1:2)], fold.raw)
        #     write.csv(sumtables.fold.raw, "Output_CellNum/SumTable_Proportions_FoldChangeRaw.csv")
        #
        #     # convert table to log2
        #
        #     fold.log2 <- log(x = fold.raw, 2)
        #
        #     sumtables.fold.log2 <- cbind(merged.df[c(1:2)], fold.log2)
        #     write.csv(sumtables.fold.log2, "Output_CellNum/SumTable_Proportions_FoldChangeLog2.csv")
        #   }
        # }
        #
        #
        # ### 4.5 - Calculate FOLD CHANGE (CELLS PER TISSUE)
        #
        # if(do.FOLDCHANGE == 1){
        #   if(use.groups == 1){
        #     if(do.CELLSPERTISSUE == 1){
        #
        #       # average of the columns for rows with specific group name
        #       sumtables.cellspertissue
        #
        #       ctrl.grp <- subset(sumtables.cellspertissue, Group == ctrls)
        #       ctrl.grp
        #
        #       ctrl.grp.means <- colMeans(ctrl.grp[-c(1:2)])
        #       as.matrix(ctrl.grp.means)
        #
        #       # As above, divide each 'cell' to the vector (by column...)
        #       fold.raw <- t(t(cellsper) / ctrl.grp.means)
        #       fold.raw
        #
        #       sumtables.fold.raw <- cbind(merged.df[c(1:2)], fold.raw)
        #       write.csv(sumtables.fold.raw, "Output_CellNum/SumTable_CellsPerTissue_FoldChangeRaw.csv")
        #
        #       # convert table to log2
        #
        #       fold.log2 <- log(x = fold.raw, 2)
        #
        #       sumtables.fold.log2 <- cbind(merged.df[c(1:2)], fold.log2)
        #       write.csv(sumtables.fold.log2, "Output_CellNum/SumTable_CellsPerTissue_FoldChangeLog2.csv")
        #     }


}







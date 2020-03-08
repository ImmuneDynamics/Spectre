##########################################################################################################
#### Spectre -- General Discovery Workflow
#### Part 2/2 - Annotation, plots, heatmaps, graphs
##########################################################################################################

    # Spectre R package: https://sydneycytometry.org.au/spectre
    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### SETUP
##########################################################################################################

    ### Load packages from library

        library(Spectre)
        Spectre::package.check() # --> change so that message at the end is "All required packages have been successfully installed"
        Spectre::package.load() # --> change so that message at the end is "All required packages have been successfully loaded"

        session_info()

        if(!require('pheatmap')) {install.packages('pheatmap')}
        if(!require('ggpubr')) {install.packages('ggpubr')}
        library(pheatmap)
        library(ggpubr)

    ### Set working directory
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

    ### Determine input directory
        setwd("Output_Spectre/Output-data/")
        InputDirectory <- getwd()
        setwd(PrimaryDirectory)

    ### Create output directory
        setwd(PrimaryDirectory)
        dir.create("Output_Spectre", showWarnings = FALSE)
        setwd("Output_Spectre")

        dir.create("Output-annotated", showWarnings = FALSE)
        setwd("Output-annotated")

        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### LOAD DATA
##########################################################################################################

    ### Read data into R
        setwd(InputDirectory)
        list.files(path = getwd(), pattern = ".csv")

        cell.dat <- fread("Clustered_TAdemo.csv")
        cell.dat

        cell.dat.sub <-fread("DimRed_TAdemo.csv")
        cell.dat.sub

    ### Read in annotation list

        setwd(PrimaryDirectory)
        setwd("metadata")

        list.files(path = getwd(), pattern = ".csv")
        annotations <- fread("annotation.csv")
        annotations

    ### Read in other stuff

        setwd(PrimaryDirectory)
        setwd("Output_Spectre/Output-info")
        list.files(path = getwd(), pattern = ".csv")

        meta.dat <- fread("metadata.csv")
        meta.dat

        CellularCols <- fread("CellularCols.csv")
        CellularCols <- as.vector(CellularCols[[1]])
        CellularCols

        ClusteringCols <- fread("ClusteringCols.csv")
        ClusteringCols <- as.vector(ClusteringCols[[1]])
        ClusteringCols

##########################################################################################################
#### ADD ANNOTATIONS
##########################################################################################################

    ### Embed population name columns in cell.dat
        cell.dat <- Spectre::do.embed.columns(x = cell.dat,
                                              type = "data.table",
                                              base.name = "FlowSOM_metacluster_42",
                                              col.name = "Population",
                                              match.to = annotations[,1],
                                              new.cols = annotations[,2])

        cell.dat

    ### Embed population name columns in cell.dat.sub

        cell.dat.sub <- Spectre::do.embed.columns(x = cell.dat.sub,
                                                 type = "data.table",
                                                 base.name = "FlowSOM_metacluster_42",
                                                 col.name = "Population",
                                                 match.to = annotations[,1],
                                                 new.cols = annotations[,2])

        cell.dat.sub

##########################################################################################################
#### GENERATE SOME NEW DATA AND SUMMARY DATA (based on population names)
##########################################################################################################

    ### Write 'large' dataset
        setwd(OutputDirectory)
        dir.create("Annotated-data", showWarnings = FALSE)
        setwd("Annotated-data")

        Spectre::write.files(dat = cell.dat,
                             file.prefix= paste0("Clustered"), # required
                             write.csv = TRUE,
                             write.fcs = TRUE)

        Spectre::write.files(dat = cell.dat,
                             file.prefix= paste0("Clustered"), # required
                             divide.by = "Sample",
                             write.csv = FALSE,
                             write.fcs = TRUE)

    ### Write 'subsample' dataset
        setwd(OutputDirectory)
        dir.create("Annotated-data", showWarnings = FALSE)
        setwd("Annotated-data")

        Spectre::write.files(dat = cell.dat.sub,
                             file.prefix = paste0("DimRed"), # required
                             write.csv = TRUE,
                             write.fcs = TRUE)

        Spectre::write.files(dat = cell.dat.sub,
                             file.prefix = paste0("DimRed"), # required
                             divide.by = "Sample",
                             write.csv = FALSE,
                             write.fcs = TRUE)

        setwd(PrimaryDirectory)


    ### Set positive cut offs for selected markers

        setwd(OutputDirectory)
        dir.create("Annotated-sumtables", showWarnings = FALSE)
        setwd("Annotated-sumtables")
        dir.create("SumTable-PercentPositive")
        setwd("SumTable-PercentPositive")

        as.matrix(names(cell.dat))

        # make.density.plot(dat = cell.dat.sub,
        #                   x.axis = "BV711.SCA.1",
        #                   y.axis = "BV605.Ly6C",
        #                   save.to.disk = TRUE)
        #
        # make.density.plot(dat = cell.dat.sub,
        #                   x.axis = "APC.BrdU",
        #                   y.axis = "BV605.Ly6C",
        #                   save.to.disk = TRUE)

        # sumtable.percent.pos(x = cell.dat,
        #                      sample.name = "SampleName",
        #                      cluster.name = "PopName",
        #                      Markers = c("BV711.SCA.1","APC.BrdU"),
        #                      Cutoff = c(580, 450)
        #                      )

    ### Generate summary data
        setwd(OutputDirectory)
        dir.create("Annotated-sumtables", showWarnings = FALSE)
        setwd("Annotated-sumtables")

        as.matrix(names(cell.dat))

        Spectre::write.sumtables(x = cell.dat,
                                 sample.col = "Sample",
                                 pop.col = "Population",
                                 measure.col = CellularCols,
                                 annot.col = names(cell.dat)[c(33:37)],
                                 group.col = "Group",
                                 do.frequencies = TRUE,
                                 cell.counts = meta.dat$`Cells per sample`,
                                 do.mfi.per.sample = TRUE,
                                 do.mfi.per.marker = TRUE)


##########################################################################################################
#### NEW PLOTS IF REQUIRED
##########################################################################################################

    ### Make a cluster plot
        setwd(OutputDirectory)
        dir.create("Annotated-plots")
        setwd("Annotated-plots")

        as.matrix(names(cell.dat.sub))

    ### Plots

        ## All data
        make.factor.plot(dat = cell.dat.sub,
                         x.axis = "UMAP_X",
                         y.axis = "UMAP_Y",
                         col.axis = "Population",
                         add.label = TRUE)


        ## Sample multi plots
        make.multi.plot(dat = cell.dat.sub,
                        x.axis = "UMAP_X",
                        y.axis = "UMAP_Y",
                        col.axis = "Population",
                        type = "factor",
                        plot.by = "Sample",
                        align.xy.by = cell.dat.sub,
                        align.col.by = cell.dat.sub)

        ## Group multi plots
        make.multi.plot(dat = cell.dat.sub,
                        x.axis = "UMAP_X",
                        y.axis = "UMAP_Y",
                        col.axis = "Population",
                        type = "factor",
                        plot.by = "Group",
                        align.xy.by = cell.dat.sub,
                        align.col.by = cell.dat.sub)


##########################################################################################################
#### PROPORTION DATA
##########################################################################################################

    ### Read in proportion and cell count data
        setwd(OutputDirectory)
        setwd("Annotated-sumtables")
        list.files(getwd())

        setwd("SumTable-Frequency/")
        list.files(getwd())

        prop <- read.csv("SumTable-Proportions.csv")

        as.matrix(names(prop))
        to.plot <- names(prop)[c(8:27)]
        prop[to.plot] <- prop[to.plot]*100
        prop

        my.comparisons <- list(c("Mock", "WNV"))

    ### Create PROPORTION AutoGraphs

        for(a in to.plot){
          # a <- "Cluster01"
          make.autograph(x = prop,
                         x.axis = "Group",
                         y.axis = a,
                         colour.by = "Batch",
                         colours = c("Black", "Red"),
                         y.axis.label = "Proportion",
                         my_comparisons = my.comparisons,
                         title = paste0("Proportion of ", a),
                         filename = paste0("Proportion of ", a, ".pdf"))
        }


    ### Create heatmaps

        prop.fold <- Spectre::do.convert.to.fold(x = prop,
                                    sample.col = "Sample",
                                    group.col = "Group",
                                    ctrl.grp = "Mock",
                                    convert.cols = c(8:27))

        prop.fold

        make.pheatmap(dat = prop.fold,
                      file.name = "Proportions.png",
                      plot.title = "Proportions",
                      sample.col = "Sample",
                      annot.cols = c(5,6),
                      rmv.cols = c(1:7),
                      is.fold = TRUE)


##########################################################################################################
#### CELL COUNT DATA
##########################################################################################################

    ### Read in proportion and cell count data
        setwd(OutputDirectory)
        setwd("Annotated-sumtables")
        list.files(getwd())

        setwd("SumTable-Frequency/")
        list.files(getwd())

        counts <- read.csv("SumTable-CellCounts.csv")

        as.matrix(names(counts))
        to.plot <- names(counts)[c(8:27)]
        counts[to.plot] <- counts[to.plot]*100
        counts

        my.comparisons <- list(c("Mock", "WNV"))

    ### Create PROPORTION AutoGraphs

        for(a in to.plot){
          # a <- "Cluster01"
          make.autograph(x = counts,
                         x.axis = "Group",
                         y.axis = a,
                         colour.by = "Batch",
                         colours = c("Black", "Red"),
                         y.axis.label = "Cells per sample",
                         my_comparisons = my.comparisons,
                         title = paste0("Number of ", a),
                         filename = paste0("Number of ", a, ".pdf"))
        }


    ### Create heatmaps

        counts.fold <- Spectre::do.convert.to.fold(x = counts,
                                                 sample.col = "Sample",
                                                 group.col = "Group",
                                                 ctrl.grp = "Mock",
                                                 convert.cols = c(8:27))

        counts.fold

        make.pheatmap(dat = counts.fold,
                      file.name = "Counts.png",
                      plot.title = "Counts",
                      sample.col = "Sample",
                      annot.cols = c(5,6),
                      rmv.cols = c(1:7),
                      is.fold = TRUE)


##########################################################################################################
#### COMING SOON - MFI HEATMAPS
##########################################################################################################

    # ### Read in proportion and cell count data
    #     setwd(OutputDirectory)
    #     setwd("Annotated-sumtables")
    #     setwd("SumTable-MFI-PerMarker")
    #     files <- list.files(getwd())
    #     files
    #
    # ### Setup
    #
    #     temp <- read.csv(files[1])
    #     as.matrix(names(temp))
    #
    #     to.plot <- c(2:13)
    #     ctrl.grp <- c(2:7)
    #
    # ### Heatmap loop
    #
    #     for(i in files){
    #       temp <- read.csv(i)
    #
    #       counts.fold <- Spectre::do.convert.to.fold(x = counts,
    #                                                  sample.col = "Sample",
    #                                                  group.col = "Group",
    #                                                  ctrl.grp = ctrl.grp,
    #                                                  convert.cols = to.plot)
    #
    #       counts.fold
    #
    #       make.pheatmap(dat = counts.fold,
    #                     file.name = "Counts.png",
    #                     plot.title = "Counts",
    #                     sample.col = "Sample",
    #                     annot.cols = c(5,6),
    #                     rmv.cols = c(1:7),
    #                     is.fold = TRUE)
    #
    #
    #
    #     }


##########################################################################################################
#### COMING SOON - PERCENT POSITIVE
##########################################################################################################























##########################################################################################################
#### PERCENT POSITIVE DATA
##########################################################################################################

### Read in and setup up FREQUENCY data

    setwd(OutputDirectory)
    setwd("SumTable-PercentPos")

    files <- list.files(path = getwd(), pattern = ".csv")
    files

    file.list <- list()

    for(i in files){
      file.list[[i]] <- read.csv(i, check.names = FALSE)
    }

    file.list[[1]]
    as.matrix(names(file.list[[1]]))

    Batches <- c(1,1,1,1,2,2,
                 1,2,2,2,1,2,
                 2,2,2,2,2) # in order of sample (air, Smoke, Smoke + FT)

    Groups <- c(rep("Air", 6),
                rep("Smoke", 6),
                rep("Smoke + FT", 5))

    for(i in names(file.list)){
      file.list[[i]][1] <- NULL
      file.list[[i]] <- cbind(file.list[[i]], Batches, Groups)
    }

    file.list[[1]]

    as.matrix(names(file.list[[1]]))
    to.plot <- c(2:25)

    ## Generic preferneces
    sample.col <- "Sample"
    group.col <- "Groups"
    ctrl.grp <- "Air"

    ## AutoGraph preferences
    file.list[[1]]
    as.matrix(names(file.list[[1]]))
    unique(file.list[[1]][[group.col]])

    to.plot <- c(2:25)
    my.comparisons <- list(c("Air", "Smoke"),
                           c("Smoke", "Smoke + FT"),
                           c("Air", "Smoke + FT"))


### Percent positive FOLD CHANGE heatmap

    setwd(OutputDirectory)
    setwd("SumTable-PercentPos")

    for(i in names(file.list)){
      # i <- names(file.list)[[1]]

      dat <- file.list[[i]]
      as.matrix(names(dat))

      dat <- convert.to.fold(x = dat,
                              sample.col = sample.col,
                              group.col = group.col,
                              ctrl.grp = ctrl.grp,
                              convert.cols = to.plot,
                              log2 = TRUE)

      i <- gsub(x = i, pattern = ".csv", "")

      as.matrix(names(dat))
      Spectre::make.pheatmap(x = dat,
                             file.name = paste0(i, ".png"),
                             plot.title = paste0(i, "Percent positive"),
                             sample.col = "Sample",
                             annot.cols = c(26,27),
                             rmv.cols = c(1),
                             row.sep = c(6,12),
                             is.fold = TRUE,
                             dendrograms = "none")
    }



### Percent positive AUTOGRAPHS

    setwd(OutputDirectory)
    setwd("SumTable-PercentPos")

    for(i in names(file.list)){
      # i <- "127I_IdU.csv"

      dat <- file.list[[i]]
      dat

      for(a in to.plot){
        # a <- 2
        a <- names(dat)[[a]]
        i <- gsub(x = i, pattern = ".csv", "")

        make.autograph(x = dat,
                       x.axis = group.col,
                       y.axis = a,
                       colour.by = "Batches",
                       colours = c("Black", "Red"),
                       y.axis.label = "Percent positive",
                       my_comparisons = my.comparisons,
                       title = paste0(i, " - % pos of ", a),
                       filename = paste0(i, " - PC pos of ", a, ".pdf")
                       )
      }
    }




##########################################################################################################
#### MFI HEATMAPS - per marker, samples vs clusters
##########################################################################################################

    setwd(OutputDirectory)
    setwd("Output_MFI_per_marker/")

    files <- list.files(path = getwd(), pattern = ".csv")
    files

    file.list <- list()

    for(i in files){
      file.list[[i]] <- read.csv(i, check.names = FALSE)
      colnames(file.list[[i]])[1] <- "Sample"
    }

    file.list[[1]]
    as.matrix(names(file.list[[1]]))

    Batches <- c(1,1,1,1,2,2,
                 1,2,2,2,1,2,
                 2,2,2,2,2) # in order of sample (air, Smoke, Smoke + FT)

    Groups <- c(rep("Air", 6),
                rep("Smoke", 6),
                rep("Smoke + FT", 5))

    for(i in names(file.list)){
      file.list[[i]] <- cbind(file.list[[i]], Batches, Groups)
    }

    file.list[[1]]

    as.matrix(names(file.list[[1]]))
    to.plot <- c(2:25)

    ## Generic preferneces
    sample.col <- "Sample"
    group.col <- "Groups"
    ctrl.grp <- "Air"

    ## AutoGraph preferences
    file.list[[1]]
    as.matrix(names(file.list[[1]]))
    unique(file.list[[1]][[group.col]])

    to.plot <- c(2:25)

### MFI heatmap loop

    setwd(OutputDirectory)
    setwd("Output_MFI_per_marker/")

    # Heatmap loop

    for(i in names(file.list)){

      dat <- file.list[[i]]

      dat <- convert.to.fold(x = dat,
                             sample.col = sample.col,
                             group.col = group.col,
                             ctrl.grp = ctrl.grp,
                             convert.cols = to.plot,
                             log2 = TRUE)


      i <- gsub(x = i, pattern = ".csv", "")

      as.matrix(names(dat))
      Spectre::make.pheatmap(x = dat,
                             file.name = paste0(i, ".png"),
                             plot.title = paste0(i, ""),
                             sample.col = c(1),
                             annot.cols = c(26,27),
                             rmv.cols = c(1),
                             row.sep = c(6,12),

                             is.fold = TRUE,
                             dendrograms = "none")
    }











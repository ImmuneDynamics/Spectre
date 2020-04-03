##########################################################################################################
#### Spectre -- General Discovery Workflow
#### Part 3/3 - Plots, heatmaps, graphs
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

    ### Determine data input directory
        InputDirectory <- PrimaryDirectory

    ### Determine metadata directory
        setwd(PrimaryDirectory)
        setwd("../../../../metadata/")
        MetaDirectory <- getwd()

    ### Create an output directory
        setwd(PrimaryDirectory)
        dir.create("Annotated-sumtables")
        setwd("Annotated-sumtables")
        OutputDirectory <- getwd()

##########################################################################################################
#### LOAD DATA
##########################################################################################################

    ### Read data into R
        setwd(InputDirectory)
        list.files(path = getwd(), pattern = ".csv")

        cell.dat <- fread("Clustered_annotated.csv")
        cell.dat

        cell.dat.sub <- fread("DimRed_annotated.csv")
        cell.dat.sub

    ### Define some parameters

        as.matrix(names(cell.dat))

        sample.name <- "Sample"
        group.name <- "Group"

        population.name <- "Population"

        to.measure <- names(cell.dat)[c(2,4:6,8:9,11:13,16:19,21:30,32)]
        to.annotate <- names(cell.dat)[c(36:37)]

##########################################################################################################
#### SUMMARY DATA
##########################################################################################################

    ### Set positive cut offs for selected markers
        setwd(OutputDirectory)
        as.matrix(names(cell.dat))

    ## Plots
        make.density.plot(dat = cell.dat.sub,
                          x.axis = "BV711.SCA.1",
                          y.axis = "BV605.Ly6C",
                          title = paste0("Density-", "SCA-1"),
                          save.to.disk = TRUE)

    ### Read in sample (cell count) metadata

        setwd(MetaDirectory)
        list.files(getwd())

        meta.dat <- read.csv("sample.details.csv")
        meta.dat

        count.vector <- as.vector(meta.dat[["Cells.per.sample"]])
        count.vector

    ### Generate summary data
        setwd(OutputDirectory)
        as.matrix(names(cell.dat))

        Spectre::write.sumtables(x = cell.dat,
                                 sample.col = sample.name,
                                 pop.col = population.name,
                                 measure.col = to.measure,
                                 annot.col = to.annotate,
                                 group.col = group.name,

                                 do.frequencies = TRUE,
                                 cell.counts = count.vector,
                                 do.mfi.per.sample = FALSE,
                                 do.mfi.per.marker = TRUE,

                                 perc.pos.markers = c("BV711.SCA.1"),
                                 perc.pos.cutoff = c(500))

##########################################################################################################
#### PROPORTION DATA
##########################################################################################################

    ### Read in proportion and cell count data

        ##
        setwd(OutputDirectory)
        setwd("SumTable-Frequency/")
        freq.list <- list.files(getwd(), ".csv")
        freq.list

        ##
        prop <- read.csv("SumTable-Proportions.csv")
        counts <- read.csv("SumTable-CellCounts.csv")

        ##
        as.matrix(names(prop))

        to.plot.nums <- c(5:24)
        to.plot <- names(prop)[to.plot.nums]
        to.annotate <- c(2:3)
        to.rmv <- c(1:4)

        prop[to.plot] <- prop[to.plot]*100
        prop

        ##
        my.comparisons <- list(c("Mock", "WNV"))
        group.colours <- c("Black", "Red")
        ctrl.group <- "Mock"

        ##
        lst <- list("Proportions" = prop, "Cells per sample" = counts)

    ### LOOP to create AutoGraphs and Pheatmaps

        for(i in names(lst)){

          dat <- lst[[i]]

          ## Autographs
          for(a in to.plot){
            # a <- "Cluster01"
            make.autograph(x = dat,
                           x.axis = group.name,
                           y.axis = a,
                           colour.by = group.name,
                           colours = group.colours,
                           y.axis.label = i,
                           my_comparisons = my.comparisons,
                           title = paste0(i, " - ", a),
                           filename = paste0(i, " - ", a, ".pdf"))
          }

          ## Pheatmaps
          dat.fold <- Spectre::do.convert.to.fold(x = dat,
                                                   sample.col = sample.name,
                                                   group.col = group.name,
                                                   ctrl.grp = ctrl.group,
                                                   convert.cols = to.plot.nums)

          dat.fold

          make.pheatmap(dat = dat.fold,
                        file.name = paste0(i, ".png"),
                        plot.title = i,
                        sample.col = sample.name,
                        annot.cols = to.annotate,
                        rmv.cols = to.rmv,
                        is.fold = TRUE)
        }



    ### Create heatmaps




##########################################################################################################
#### CELL COUNT DATA
##########################################################################################################

    ### Read in proportion and cell count data
        setwd(InputDirectory)
        setwd("SumTable-Frequency/")
        list.files(getwd(), ".csv")

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
#### MFI (per marker)
##########################################################################################################

    ### Read in proportion and cell count data
        setwd(InputDirectory)
        setwd("SumTable-MFI-PerMarker")
        files <- list.files(getwd(), ".csv")
        files

    ### Setup

        setup <- read.csv(files[1])
        as.matrix(names(setup))

        to.plot <- c(8:27)
        to.annot <- c(6:7)
        to.rmv <- c(1:7)

        samp.col <- "Sample"
        grp.col <- "Group"
        ctrl.grp <- "Mock"

        my.comparisons <- list(c("Mock", "WNV"))

    ### Create  AutoGraphs

        plot.names <- names(setup)[to.plot]

        for(i in files){
          dat <- read.csv(i)

          nme <- gsub(".csv", "", i)
          nme <- gsub("SumTable-MFI-*", "", nme)

          for(a in plot.names){
            make.autograph(x = dat,
                           x.axis = grp.col,
                           y.axis = a,
                           colour.by = "Batch",
                           colours = c("Black", "Red"),
                           y.axis.label = "MFI",
                           my_comparisons = my.comparisons,
                           title = paste0(nme, " MFI", " for ", a),
                           filename = paste0(nme, " MFI", " for ", a, ".pdf"))
          }
        }

    ### Heatmap loop

        for(i in files){
          temp <- read.csv(i)

          nme <- gsub(".csv", "", i)
          nme <- gsub("SumTable-MFI-*", "", nme)

          counts.fold <- Spectre::do.convert.to.fold(x = temp,
                                                     sample.col = samp.col,
                                                     group.col = grp.col,
                                                     ctrl.grp = ctrl.grp,
                                                     convert.cols = to.plot)

          counts.fold

          make.pheatmap(dat = counts.fold,
                        file.name = paste0(nme, ".png"),
                        plot.title = paste0(nme),
                        sample.col = samp.col,
                        annot.cols = to.annot,
                        rmv.cols = to.rmv,
                        fold.max.range = 2,
                        fold.min.range = -2,
                        row.sep = 6,
                        dendrograms = "none",
                        is.fold = TRUE)
        }


##########################################################################################################
#### PERCENT POSITIVE
##########################################################################################################

    ### Read in proportion and cell count data
        setwd(InputDirectory)
        setwd("SumTable-MFI-PercentPositive")
        files <- list.files(getwd(), ".csv")
        files

    ### Setup

        setup <- read.csv(files[1], check.names = FALSE)
        as.matrix(names(setup))

        to.plot <- c(8:27)
        to.annot <- c(5:6)
        to.rmv <- c(1:7)

        samp.col <- "Sample"
        grp.col <- "Group"
        ctrl.grp <- "Mock"

        my.comparisons <- list(c("Mock", "WNV"))

    ### Create  AutoGraphs

        plot.names <- names(setup)[to.plot]

        for(i in files){
          dat <- read.csv(i)

          nme <- gsub(".csv", "", i)
          nme <- gsub("SumTable-PercPos-*", "", nme)

          for(a in plot.names){
            make.autograph(x = dat,
                           x.axis = grp.col,
                           y.axis = a,
                           colour.by = "Batch",
                           colours = c("Black", "Red"),
                           y.axis.label = "Percent positive",
                           my_comparisons = my.comparisons,
                           title = paste0(a, " - ", nme, " percent positive"),
                           filename = paste0(nme, " - percent positive of ", a, ".pdf"))
          }
        }

    ### Heatmap loop

        for(i in files){
          temp <- read.csv(i)

          nme <- gsub(".csv", "", i)
          nme <- gsub("SumTable-PercPos-*", "", nme)

          counts.fold <- Spectre::do.convert.to.fold(x = temp,
                                                     sample.col = samp.col,
                                                     group.col = grp.col,
                                                     ctrl.grp = ctrl.grp,
                                                     convert.cols = to.plot)

          counts.fold

          make.pheatmap(dat = counts.fold,
                        file.name = paste0(nme, ".png"),
                        plot.title = paste0(nme),
                        sample.col = samp.col,
                        annot.cols = to.annot,
                        rmv.cols = to.rmv,
                        fold.max.range = 2,
                        fold.min.range = -2,
                        row.sep = 6,
                        dendrograms = "none",
                        is.fold = TRUE)
        }



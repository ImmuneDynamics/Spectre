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

    ### Determine input directory
        setwd("Output_Spectre/Annotated-sumtables")
        InputDirectory <- getwd()

##########################################################################################################
#### PROPORTION DATA
##########################################################################################################

    ### Read in proportion and cell count data
        setwd(InputDirectory)
        setwd("SumTable-Frequency/")
        list.files(getwd(), ".csv")

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
#### COMING SOON - MFI (per marker) HEATMAPS
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
#### COMING SOON - PERCENT POSITIVE
##########################################################################################################


    ### Read in proportion and cell count data
        setwd(InputDirectory)
        setwd("SumTable-PercentPositive")
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






















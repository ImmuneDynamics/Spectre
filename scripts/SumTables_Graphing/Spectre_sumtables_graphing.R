##########################################################################################################
#### Spectre -- SumTables, Heatmaphs, Graphs
##########################################################################################################

    # Spectre R package: https://sydneycytometry.org.au/spectre
    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### SETUP
##########################################################################################################

    ### Load packages from library

        library(Spectre)
        Spectre::package.check()
        Spectre::package.load()

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

    ### Create an output directory
        setwd(PrimaryDirectory)
        dir.create("Output-sumtables")
        setwd("Output-sumtables")
        OutputDirectory <- getwd()

##########################################################################################################
#### LOAD DATA
##########################################################################################################

    ### Read data into R
        setwd(InputDirectory)
        list.files(path = getwd(), pattern = ".csv")

        cell.dat <- fread("Clustered_TAdemo.csv")
        cell.dat

    ### Read in any metadata
        setwd(InputDirectory)
        setwd("metadata/")

        list.files(path = getwd(), pattern = ".csv")

        meta.dat <- fread("sample.details.csv")
        meta.dat

    ### Choose columns to measure

        as.matrix(names(cell.dat))

        to.measure <- names(cell.dat)[c(2,4:6,8:9,11:13,16:19,21:30,32)]
        to.measure

    ### Define some parameters

        as.matrix(names(cell.dat))

        sample.col <- "Sample"
        group.col <- "Group"
        batch.col <- "Batch"
        pop.col <- "FlowSOM_metacluster"

    ### Define groups and statistical comparisons

        as.matrix(unique(cell.dat[[group.col]]))

        ctrl.grp <- "Mock"

        grp.order <- c("Mock", "WNV")
        grp.colours <- c("Black", "Red")

        # TODO remove no longer required as autograph is not supporting it
        # stat.comparisons <- list(c("Mock", "WNV")) # A list of comparisons for statistical test (used in graphing and stats)
        # 
        # var.test <- "kruskal.test" # can be "kruskal.test", "anova", or NULL
        # pair.test <- "wilcox.test" # can be "wilcox.test". "t.test", or NULL

    ### Define cell counts (if desired)

        meta.dat
        as.matrix(unique(cell.dat[[sample.col]]))

        cell.counts <- as.vector(meta.dat[["Cells per sample"]])
        cell.counts

##########################################################################################################
#### Plots to aid in summary data generation
##########################################################################################################

    ### Set positive cut offs for selected markers

        as.matrix(names(cell.dat))

        plot.dat <- do.subsample(cell.dat,
                                 method = "per.sample",
                                 samp.col = group.col,
                                 targets = rep(5000, length(grp.order)))

        plot.y <- "BV605.Ly6C"
        to.plot <- c("BV711.SCA.1", "APC.BrdU")

    ### Create some plots

        setwd(OutputDirectory)

        for(i in to.plot){
          make.multi.plot(dat = plot.dat,
                         x.axis = i,
                         y.axis = plot.y,
                         col.axis = group.col,
                         type = "factor",
                         plot.by = group.col,
                         figure.title = paste0(i, " - split by group"),
                         save.each.plot = TRUE)
        }

    ### Define cutoffs
        as.matrix(names(cell.dat))

        markers.cutoff <- c("BV711.SCA.1", "APC.BrdU")
        values.cutoff <- c(600, 400)

#########################################################################################################
#### Create summary data (per non-annotated cluster) and produce graphs and heatmaps
#########################################################################################################

    ### Write sumtables - proportions, cell counts, MFI
        setwd(OutputDirectory)

        write.sumtables(dat = cell.dat,
                        sample.col = sample.col,
                        pop.col = pop.col,

                        measure.col = to.measure,
                        annot.col = c(group.col, batch.col),
                        group.col = group.col,

                        cell.counts = cell.counts, # vector must be in order of the samples in which they appear (unique(cell.dat[[sample.col]]))
                        do.mfi.per.sample = FALSE,
                        do.mfi.per.marker = TRUE)


    ### Write sumtables for 'percent positive' only
        setwd(OutputDirectory)

        write.sumtables(dat = cell.dat,
                        sample.col = sample.col,
                        pop.col = pop.col,

                        measure.col = to.measure,
                        annot.col = c(group.col, batch.col),
                        group.col = group.col,
                        
                        do.mfi.per.sample = FALSE,
                        do.mfi.per.marker = FALSE,

                        perc.pos.markers = markers.cutoff,
                        perc.pos.cutoff = values.cutoff)

#########################################################################################################
#### Produce graphs and heatmaps
#########################################################################################################

    ### List of sumtables
        setwd(OutputDirectory)

        sumtable.files <- list.files(getwd(), ".csv")
        sumtable.files

    ### Select column names to plot from sumtable files

        temp <- fread(sumtable.files[[1]])
        as.matrix(names(temp))

        plot.names <- names(temp)[c(6:20)]
        plot.names <- as.character(plot.names)

    ### Pheatmap loop
        setwd(OutputDirectory)

        for(i in sumtable.files){
          dat <- fread(i)

          ## Convert to fold
          dat.fold <- do.convert.to.fold(x = dat,
                                         sample.col = sample.col,
                                         group.col = group.col,
                                         ctrl.grp = ctrl.grp,
                                         convert.cols = plot.names)

          ## Remove "Inf" or -Inf"
          dat.fold[sapply(dat.fold, is.infinite)] <- NA

          ## Make Pheatmap
          a <- gsub(".csv", "", i)

          make.pheatmap(dat = dat.fold,
                        file.name = paste0(a, ".png"),
                        plot.title = a,
                        sample.col = sample.col,
                        annot.cols = group.col,
                        plot.cols = plot.names,
                        dendrograms = "none",
                        is.fold = TRUE)
        }


    ### AutoGraph loops
        setwd(OutputDirectory)

        for(i in sumtable.files){
          dat <- fread(i)

          if(grepl('Cells per', dat[1,1], fixed = TRUE)){
            scale <- "sci"
          }

          if(!grepl('Cells per', dat[1,1], fixed = TRUE)){
            scale <- "lin"
          }

          for(a in plot.names){
            make.autograph(dat = dat,
                           x.axis = group.col,
                           grp.order = grp.order,
                           y.axis = a,
                           colour.by = group.col,
                           colours = grp.colours,
                           y.axis.label = dat[1,1],
                           title = paste0(a),
                           scale = scale,
                           filename = paste0(dat[1,1], " - ", a, ".pdf"))
          }
        }


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

        cell.dat.sub <- fread("DimRed_TAdemo.csv")
        cell.dat.sub

    ### Read in any metadata
        setwd(InputDirectory)
        list.files(path = getwd(), pattern = ".csv")

        meta.dat <- fread("sample.details.csv")
        meta.dat

    ### Define some parameters

        sample.col <- "Sample"
        group.col <- "Group"
        batch.col <- "Batch"
        pop.col <- "FlowSOM_metacluster"

        to.measure <- names(cell.dat)[c(2,4:6,8:9,11:13,16:19,21:30,32)]

        ctrl.grp <- "Mock"
        my.comparisons <- list(c("Mock", "WNV")) # A list of comparisons for statistical test (used in graphing and stats)

        cell.counts <- as.vector(meta.dat[["Cells per sample"]])
        cell.counts

##########################################################################################################
#### SUMMARY DATA
##########################################################################################################

    ### Set positive cut offs for selected markers
        setwd(OutputDirectory)
        as.matrix(names(cell.dat))

    ## Plots
        make.factor.plot(dat = cell.dat.sub,
                          x.axis = "BV711.SCA.1",
                          y.axis = "BV605.Ly6C",
                          col.axis = group.col,
                          title = paste0("SCA-1"),
                          save.to.disk = TRUE)

#########################################################################################################
#### Create summary data (per non-annotated cluster) and produce graphs and heatmaps
#########################################################################################################

    ### Write sumtables
        setwd(OutputDirectory)

        write.sumtables(x = cell.dat,
                        sample.col = sample.col,
                        pop.col = pop.col,

                        measure.col = to.measure,
                        annot.col = c(group.col, batch.col),
                        group.col = group.col,

                        do.frequencies = TRUE,
                        cell.counts = cell.counts,
                        do.mfi.per.sample = FALSE,
                        do.mfi.per.marker = TRUE)

#########################################################################################################
#### Produce graphs and heatmaps
#########################################################################################################

    ### List of sumtables
        setwd(OutputDirectory)

        sumtable.files <- list.files(getwd(), ".csv")
        sumtable.files

    ### Plot names

        plot.names <- unique(cell.dat[[pop.col]])
        plot.names <- sort(plot.names, decreasing = FALSE)
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
                                         convert.cols = plot.names) ########

          a <- gsub(".csv", "", i)

          ## Make Pheatmap
          make.pheatmap(dat = dat.fold,
                        file.name = paste0(a, ".png"),
                        plot.title = a,
                        sample.col = sample.col,
                        annot.cols = group.col,
                        plot.cols = plot.names,
                        fold.range = c(2, -2),
                        dendrograms = "none",
                        is.fold = TRUE)
        }


    ### AutoGraph loops
        setwd(OutputDirectory)

        for(i in sumtable.files){
          dat <- fread(i)

          for(a in plot.names){
            make.autograph(x = dat,
                           x.axis = group.col,
                           y.axis = a,
                           colour.by = group.col,
                           colours = c("Black", "Red"),
                           y.axis.label = dat[1,1],
                           my_comparisons = my.comparisons,
                           title = paste0(a),
                           filename = paste0(dat[1,1], " - ", a, ".pdf"))
          }
        }



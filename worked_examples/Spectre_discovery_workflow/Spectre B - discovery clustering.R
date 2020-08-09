##########################################################################################################
#### Spectre discovery workflow - B - discovery clustering
##########################################################################################################

    # Spectre R package: https://sydneycytometry.org.au/spectre
    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### Analysis session setup
##########################################################################################################

    ### Load packages

        library(Spectre)
        Spectre::package.check()    # Check that all required packages are installed
        Spectre::package.load()     # Load required packages

    ### Set DT threads

        getDTthreads()

    ### Set primary directory

        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

    ### Set input directory

        setwd(PrimaryDirectory)
        setwd("Output A - data prep/")
        InputDirectory <- getwd()
        InputDirectory
        setwd(PrimaryDirectory)

    ### Set metadata directory

        setwd(PrimaryDirectory)
        setwd("metadata/")
        MetaDirectory <- getwd()
        MetaDirectory
        setwd(PrimaryDirectory)

    ### Set output directory

        setwd(PrimaryDirectory)
        dir.create("Output B - discovery analysis", showWarnings = FALSE)
        setwd("Output B - discovery analysis")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### Import data
##########################################################################################################

    setwd(InputDirectory)

    ### Import data

        cell.dat <- fread('Cellular.data.csv')
        cell.dat

    ### Import preferences

        prefs <- readRDS(file = 'Analysis preferences.rds')
        prefs

##########################################################################################################
#### Run clustering and DR
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output - clustered")
    setwd("Output - clustered")

    ### Run clustering

        cell.dat <- run.flowsom(cell.dat, prefs$clustering.cols, meta.k = prefs$metak)
        cell.dat

        fwrite(cell.dat, "Clustered.csv")

    ### Run DR

        cell.sub <- do.subsample(cell.dat, prefs$sub.targets, prefs$sub.by)
        cell.sub <- run.umap(cell.sub, prefs$clustering.cols)
        cell.sub

        fwrite(cell.dat, "Clustered.RD.csv")

    ### Make expression heatmap

        exp <- do.aggregate(cell.dat, prefs$cellular.cols, by = "FlowSOM_metacluster")
        exp

        make.pheatmap(exp, "FlowSOM_metacluster", plot.cols = prefs$cellular.cols)

    ### Make expression plots

        make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE)
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", prefs$group.col, col.type = 'factor')
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", prefs$cellular.cols)

        for(i in prefs$cellular.cols){
          make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", i, group.col, figure.title = paste0('Multiplot - ', i, '.png'))
        }

##########################################################################################################
#### Annotate clusters and write summary data
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output - annotated")
    setwd("Output - annotated")

    ### Identify cellular populations

        annots <- list("CD4 T cells" = c(6),
                       "CD8 T cells" = c(5),
                       "NK cells" = c(4),
                       "Neutrophils" = c(1),
                       "Infil Macrophages" = c(3),
                       "Microglia" = c(2,7)
                       )

        annots <- do.list.switch(annots)
        setorderv(annots, 'Values')
        annots

    ### Add population names to datasets

        names(annots) <- c('Values', "Population")
        annots

        cell.dat <- do.add.cols(cell.dat, "FlowSOM_metacluster", annots, "Values")
        cell.dat

        cell.sub <- do.add.cols(cell.sub, "FlowSOM_metacluster", annots, "Values")
        cell.sub

    ### Save annotated data and make population plots

        fwrite(cell.dat, "Clustered.annotated.csv")
        fwrite(cell.sub, "Clustered.DR.annotated.csv")

    ### Population plots

        annot.exp <- do.aggregate(cell.dat, prefs$cellular.cols, by = "Population")
        make.pheatmap(annot.exp, "Population", plot.cols = prefs$cellular.cols)

        make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "Population", col.type = 'factor', add.label = TRUE)
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "Population", prefs$group.col, col.type = 'factor')

    ### Columns for differential expression

        setwd(OutputDirectory)
        dir.create("Output - annotated")
        setwd("Output - annotated")

        dir.create("Percent positive plots")
        setwd("Percent positive plots")

        as.matrix(names(cell.dat))

        prefs$perc.pos.markers <- names(cell.dat)[c(18)]

        for(i in prefs$perc.pos.markers){
          make.multi.plot(cell.sub, i, "CD45_asinh", plot.by = prefs$group.col)
        }

        prefs$perc.pos.cutoffs <- c(1.75)

        prefs$perc.pos.markers
        prefs$perc.pos.cutoffs

        for(i in c(1:length(prefs$perc.pos.markers))){
          a <- prefs$perc.pos.markers[[i]]
          b <- prefs$perc.pos.cutoffs[[i]]

          x <- cell.sub[cell.sub[[a]] >= b]
          x$Pos <- TRUE
          y <- cell.sub[cell.sub[[a]] < b]
          y$Pos <- FALSE
          z <- rbind(x,y)

          make.multi.plot(z, a, "CD45_asinh", 'Pos', prefs$group.col)
        }

##########################################################################################################
#### Write summary data
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output - summary data")
    setwd("Output - summary data")

    ### Setup cell count data

        as.matrix(unique(cell.dat[[prefs$sample.col]]))
        counts <- prefs$sample.dat$`Cells per sample`

    ### Create summary tables

        write.sumtables(cell.dat,
                        sample.col = prefs$sample.col,
                        group.col = prefs$group.col,
                        pop.col = "Population",
                        measure.col = prefs$cellular.cols,
                        annot.col = c(prefs$group.col, prefs$batch.col),
                        do.proportions = TRUE,
                        cell.counts = counts,
                        do.mfi.per.sample = FALSE,
                        do.mfi.per.marker = TRUE,
                        perc.pos.markers = prefs$perc.pos.markers,
                        perc.pos.cutoff = prefs$perc.pos.cutoffs)

##########################################################################################################
#### Spectre: General Discovery Workflow - (3/4) - Clustering and dimensionality reduction
##########################################################################################################

    # Spectre R package: https://github.com/ImmuneDynamics/Spectre
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
        setwd("Output 2 - batch alignment/Output 4 - fine alignment/F - Fine aligned data/")
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
        dir.create("Output 3 - clustering and DR", showWarnings = FALSE)
        setwd("Output 3 - clustering and DR")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### Import data
##########################################################################################################

    setwd(InputDirectory)

    ### Import data

        list.files(getwd(), '.csv')
        
        cell.dat <- fread('cell.dat.csv')
        cell.dat

##########################################################################################################
#### Import meta data
##########################################################################################################

    setwd(MetaDirectory)
        
    ### Import metadata
        
        list.files(getwd(), '.csv')
        
        meta.dat <- fread('sample.details.csv')
        meta.dat
    
##########################################################################################################
#### Setup preferences
##########################################################################################################
        
    ### Sample preferences

        as.matrix(names(cell.dat))

        sample.col <- "Sample"
        group.col <- "Group"
        batch.col <- "Batch"

    ### Downsample preferences
        
        sub.by <- "Group"
        
        as.matrix(unique(cell.dat[[sub.by]]))

        cells.per <- list()
        
        for(i in unique(cell.dat[[sub.by]])){
            # i <- unique(cell.dat[[sub.by]])[[1]]
            cells.per[[i]] <- nrow(cell.dat[cell.dat[[sub.by]] == i,])
        }
        
        cells.per
        
        sub.targets <- c(2000, 18380)
        sub.targets
        
    ### Clustering preferences

        ## Cellular cols
        as.matrix(names(cell.dat))

        cellular.cols <- names(cell.dat)[c(33:41)]
        cellular.cols

        ## Columns for clustering
        as.matrix(names(cell.dat))

        clustering.cols <- names(cell.dat)[c(33:41)]
        clustering.cols

        ## Cluster numbers etc
        metak <- 'auto'

##########################################################################################################
#### Run clustering and dimensionality reduction
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 1 - clustered")
    setwd("Output 1 - clustered")

    ### Run clustering

        cell.dat <- run.flowsom(cell.dat, clustering.cols, meta.k = metak)
        cell.dat

        fwrite(cell.dat, "Clustered.csv")

    ### Run DR

        cell.sub <- do.subsample(cell.dat, sub.targets, sub.by)
        cell.sub <- run.umap(cell.sub, clustering.cols)
        cell.sub

        fwrite(cell.sub, "RD.sub.csv")

    ### Make expression heatmap

        exp <- do.aggregate(cell.dat, cellular.cols, by = "FlowSOM_metacluster")
        exp

        make.pheatmap(exp, "FlowSOM_metacluster", plot.cols = cellular.cols)

    ### Make expression plots

        make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE)
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", group.col, col.type = 'factor')
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", cellular.cols)

        for(i in cellular.cols){
          make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", i, group.col, figure.title = paste0('Multiplot - ', i, '.png'))
        }

##########################################################################################################
#### Annotate clusters and write summary data
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 2 - annotated")
    setwd("Output 2 - annotated")

    ### Identify cellular populations

        annots <- list("CD4 T cells" = c(7),
                       "CD8 T cells" = c(6),
                       "NK cells" = c(5),
                       "Neutrophils" = c(1),
                       "Infil Macrophages" = c(3),
                       "Microglia" = c(2,4,8)
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
        fwrite(cell.sub, "DR.sub.annotated.csv")

    ### Population plots

        annot.exp <- do.aggregate(cell.dat, cellular.cols, by = "Population")
        make.pheatmap(annot.exp, "Population", plot.cols = cellular.cols)

        make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "Population", col.type = 'factor', add.label = TRUE)
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "Population", group.col, col.type = 'factor')

##########################################################################################################
#### Percent positive plots and cutoffs
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 2 - annotated")
    setwd("Output 2 - annotated")
    
    dir.create("Percent positive plots")
    setwd("Percent positive plots")
        
    ### Plots for cytoff

        as.matrix(names(cell.dat))

        perc.pos.markers <- names(cell.dat)[c(40)]
        plot.against <- 'CD11b_asinh_fineAlign'
        
        for(i in perc.pos.markers){
          make.multi.plot(cell.sub, i, plot.against, plot.by = group.col)
        }

    ### Percent positive cutoffs
        
        perc.pos.markers
        perc.pos.cutoffs <- c(1.75)
        
        perc.pos.markers
        perc.pos.cutoffs
        
        for(i in c(1:length(perc.pos.markers))){
          a <- perc.pos.markers[[i]]
          b <- perc.pos.cutoffs[[i]]

          x <- cell.sub[cell.sub[[a]] >= b]
          x$Pos <- TRUE
          y <- cell.sub[cell.sub[[a]] < b]
          y$Pos <- FALSE
          z <- rbind(x,y)

          make.multi.plot(z, a, plot.against, 'Pos', group.col)
        }

    ### Generate data.table
        
        perc.pos <- data.table(perc.pos.markers, perc.pos.cutoffs)
        names(perc.pos) <- c('Columns', 'Cutoff')
        perc.pos
                
##########################################################################################################
#### Write summary data
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 3 - summary data")
    setwd("Output 3 - summary data")
    
    ### Select columns to measure MFI
    
        as.matrix(cellular.cols)
        dyn.cols <- cellular.cols[c(5,8)]
        dyn.cols
    
    ### Setup cell count data

        as.matrix(unique(cell.dat[[sample.col]]))
        meta.dat
        
        counts <- meta.dat[,c(sample.col, 'Cells per sample'), with = FALSE]
        counts

    ### Create summary tables

        sum.dat <- create.sumtable(dat = cell.dat, 
                                   sample.col = sample.col,
                                   pop.col = "Population",
                                   use.cols = dyn.cols, 
                                   annot.cols = c(group.col, batch.col), 
                                   counts = counts, 
                                   perc.pos = perc.pos)
    
        as.matrix(names(sum.dat))
        sum.dat    

    ### Write summary data
        
        fwrite(sum.dat, 'sum.dat.csv')
    
##########################################################################################################
#### Output session info
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output - info")
    setwd("Output - info")
    
    sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
    session_info()
    sink()
    
    #saveRDS(prefs, "Analysis preferences.rds") # analysis preferences
        
        
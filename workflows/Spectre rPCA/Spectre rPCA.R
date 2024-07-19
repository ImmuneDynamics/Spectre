##########################################################################################################
#### Spectre - rPCA Batch Integration Workflow
##########################################################################################################

    # Spectre R package: https://immunedynamics.io/spectre
    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### Create a folder structure for your analysis run
##########################################################################################################

    ### Create a master folder with a meaningful name. Then inside that folder, insert the following:
    
        # One folder called 'data' -- this will contain your data CSV or FCS files
        # One folder called 'metadata' -- this will contain a CSV containg your sample metadata
        # One folder called 'Spectre rPCA' or similar -- place this analysis script there
    
    ### Example:
    
        # BM analysis
        #   /data
        #       -- Contains data files, one CSV or FCS per sample
        #   /metadata
        #       -- Contains a CSV containing sample metadata (group, batch, etc)
        #   /Spectre rPCA
        #       -- Spectre rPCA.R

##########################################################################################################
#### 1. Load packages, and set working directory
##########################################################################################################

    ### Load libraries

        library(Spectre)
        Spectre::package.check()    # Check that all required packages are installed
        Spectre::package.load()     # Load required packages

    ### Install Seurat package
        
        if(!require('Seurat')) {install.packages('Seurat')}
        library('Seurat')
        
    ### Set PrimaryDirectory
        
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory
        
    ### Set 'input' directory
        
        setwd(PrimaryDirectory)
        dir.create('../data', showWarnings = FALSE)
        setwd("../data/")
        InputDirectory <- getwd()
        setwd(PrimaryDirectory)

    ### Set 'metadata' directory
        
        setwd(PrimaryDirectory)
        dir.create('../metadata', showWarnings = FALSE)
        setwd("../metadata/")
        MetaDirectory <- getwd()
        setwd(PrimaryDirectory)

    ### Create output directory
        
        dir.create("Output_Spectre", showWarnings = FALSE)
        setwd("Output_Spectre")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### Demo dataset
##########################################################################################################

    ### If you need the demo dataset, uncomment the following code (select all, CMD+SHIFT+C) and run to download
    ### Alternative: download from https://github.com/ImmuneDynamics/data/blob/main/simBatches.zip?raw=TRUE

        # setwd(PrimaryDirectory)
        # setwd("../")
        # getwd()
        # download.file(url = "https://github.com/ImmuneDynamics/data/blob/main/simBatches.zip?raw=TRUE", destfile = 'simBatches.zip', mode = 'wb')
        # unzip(zipfile = 'simBatches.zip')
        # for(i in list.files('simBatches/data', full.names = TRUE)){
        #   file.rename(from = i,  to = gsub('simBatches/', '', i))
        # }
        # for(i in list.files('simBatches/metadata', full.names = TRUE)){
        #   file.rename(from = i,  to = gsub('simBatches/', '', i))
        # }
        # unlink(c('simBatches/', 'simBatches.zip', '__MACOSX'), recursive = TRUE)
        
##########################################################################################################
#### 2. Import and prep data
##########################################################################################################

    ### Import data

        setwd(InputDirectory)
        list.files(InputDirectory, ".csv")

        data.list <- Spectre::read.files(file.loc = InputDirectory,
                                         file.type = ".csv",
                                         do.embed.file.names = TRUE)

    ### Check the data

        check <- do.list.summary(data.list)

        check$name.table # Review column names and their subsequent values
        check$ncol.check # Review number of columns (features, markers) in each sample
        check$nrow.check # Review number of rows (cells) in each sample

        data.list[[1]]

    ### Merge data

        cell.dat <- Spectre::do.merge.files(dat = data.list)
        cell.dat

    ### Read in metadata  
       
        setwd(MetaDirectory)
        
        meta.dat <- fread("sample.details.csv")
        meta.dat
        
##########################################################################################################
#### 3. Data transformation
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 1 - transformed plots")
    setwd("Output 1 - transformed plots")
        
    ### Arcsinh transformation

        as.matrix(names(cell.dat))

        to.asinh <- names(cell.dat)[c(1:8)]
        to.asinh

        cofactor <- 500
        plot.against <- "BV605 Ly6C_asinh"

        cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
        transformed.cols <- paste0(to.asinh, "_asinh")

        for(i in transformed.cols){
          make.colour.plot(do.subsample(cell.dat, 20000), i, plot.against)
        }

##########################################################################################################
#### 4. Add metadata and set some preferences
##########################################################################################################

    ### Add metadata to data.table

        meta.dat
        sample.info <- meta.dat[,c(1:4)]
        sample.info
        
        meta.dat
        
        cell.dat <- do.add.cols(cell.dat, "FileName", sample.info, "FileName", rmv.ext = TRUE)
        cell.dat

    ### Cellular columns

        as.matrix(names(cell.dat))

        cellular.cols <- names(cell.dat)[c(11:18)]
        as.matrix(cellular.cols)

    ### Clustering columns        
        
        as.matrix(names(cell.dat))
        
        cluster.cols <- names(cell.dat)[c(11:18)]
        as.matrix(cluster.cols)

    ### Clusters to align
        
        as.matrix(names(cell.dat))
        
        raw.cols <- names(cell.dat)[c(11:18)]
        as.matrix(raw.cols)

    ### Factors
        
        exp.name <- "BM experiment"
    
        as.matrix(names(cell.dat))
        
        sample.col <- "Sample"
        group.col <- "Group"
        batch.col <- "Batch"

    ### Subsample targets per group

        data.frame(table(cell.dat[[group.col]])) # Check number of cells per sample.

        unique(cell.dat[[group.col]])
        
        sub.targets.group <- c(20000, 20000) # target subsample numbers from each group
        sub.targets.group
        
    ### Subsample targets per batch
        
        data.frame(table(cell.dat[[batch.col]])) # Check number of cells per sample.
        
        unique(cell.dat[[batch.col]])
        
        sub.targets.batch <- c(20000, 20000) # target subsample numbers from each group
        sub.targets.batch

    ### Choose a batch as reference
        
        as.matrix(unique(cell.dat[[batch.col]]))
        
        ref <- 'A'
        ref
        
##########################################################################################################
#### 5. Testing batch integration
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 2 - Test batch integration")
    setwd("Output 2 - Test batch integration")

    ### Pre-alignment assessment
    
        setwd(OutputDirectory)
        setwd("Output 2 - Test batch integration")
        dir.create("2.1 - Pre-alignment")
        setwd("2.1 - Pre-alignment")
        
        as.matrix(unique(cell.dat[[batch.col]]))
        
        test <- do.subsample(cell.dat, sub.targets.batch, divide.by = batch.col)
        test
        
        test <- run.umap(test, raw.cols)
        test
        
        make.colour.plot(test, 'UMAP_X', 'UMAP_Y', batch.col)
        make.multi.plot(test, 'UMAP_X', 'UMAP_Y', cellular.cols)
        make.multi.plot(test, 'UMAP_X', 'UMAP_Y', batch.col, batch.col)
    
    ### Alignment test
        
        setwd(OutputDirectory)
        setwd("Output 2 - Test batch integration")
        dir.create("2.2 - Test alignment")
        setwd("2.2 - Test alignment")
        
        test <- run.rpca(dat = test, use.cols = raw.cols, batch.col = batch.col, reference = ref)
        test

        aligned.cols <- paste0(raw.cols, '_rPCA_aligned')
        aligned.cols
        
        test[,..aligned.cols]
        
        test <- run.umap(test, aligned.cols, umap.x.name = 'UMAP_X_Integrated', umap.x.name = 'UMAP_Y_Integrated')
        test
        
        make.colour.plot(test, 'UMAP_X_Integrated', 'UMAP_Y_Integrated', batch.col)
        
        make.multi.plot(test, 'UMAP_X_Integrated', 'UMAP_Y_Integrated', raw.cols, figure.title = 'Raw cols')
        make.multi.plot(test, 'UMAP_X_Integrated', 'UMAP_Y_Integrated', aligned.cols, figure.title = 'Aligned cols')
        make.multi.plot(test, 'UMAP_X_Integrated', 'UMAP_Y_Integrated', batch.col, batch.col)
        
    ### Comparison plots

        for(i in raw.cols){
          make.colour.plot(test, paste0(i, '_rPCA_aligned'), i, batch.col)
        }
        
    ### Clean up
        
        rm(aligned.cols)
        
##########################################################################################################
#### 6. Full batch integration
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 3 - Full batch integration")
    setwd("Output 3 - Full batch integration")
    
    ### Batch integration
        
        cell.dat <- run.rpca(dat = cell.dat, use.cols = raw.cols, batch.col = batch.col, reference = ref)
        cell.dat
            
        aligned.cols <- paste0(raw.cols, '_rPCA_aligned')
        aligned.cols

        cell.dat[,..aligned.cols]
        
        batch.sub <- do.subsample(cell.dat, sub.targets.batch, divide.by = batch.col)
        batch.sub
        
        batch.sub <- run.umap(batch.sub, aligned.cols, UMAP.x.name = 'UMAP_X_Integrated', UMAP.y.name = 'UMAP_Y_Integrated')
        batch.sub
        
    ### Plots
        
        make.colour.plot(batch.sub, 'UMAP_X_Integrated', 'UMAP_Y_Integrated', batch.col)
        
        make.multi.plot(batch.sub, 'UMAP_X_Integrated', 'UMAP_Y_Integrated', raw.cols, figure.title = 'Raw cols')
        make.multi.plot(batch.sub, 'UMAP_X_Integrated', 'UMAP_Y_Integrated', aligned.cols, figure.title = 'Aligned cols')
        make.multi.plot(batch.sub, 'UMAP_X_Integrated', 'UMAP_Y_Integrated', batch.col, batch.col)
    
    ### Comparison plots
        
        for(i in raw.cols){
          make.colour.plot(batch.sub, paste0(i, '_rPCA_aligned'), i, batch.col)
        }
        
##########################################################################################################
#### 7. Clustering and dimensionality reduction
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 4 - clustering")
    setwd("Output 4 - clustering")

    ### Update clustering columns
    
        cellular.cols <- paste0(cellular.cols, '_rPCA_aligned')
        cellular.cols
    
        cluster.cols <- paste0(cluster.cols, '_rPCA_aligned')
        cluster.cols
    
    ### Clustering

        cell.dat <- run.flowsom(cell.dat, cluster.cols, meta.k = 15)
        cell.dat
        
    ### Dimensionality reduction

        cell.sub <- do.subsample(cell.dat, sub.targets.group, group.col)
        cell.sub
        
        cell.sub <- run.umap(cell.sub, cluster.cols)
        cell.sub

    ### DR plots

        make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE)
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", cellular.cols)
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", group.col, col.type = 'factor')

    ### Expression heatmap

        exp <- do.aggregate(cell.dat, cellular.cols, by = "FlowSOM_metacluster")
        make.pheatmap(exp, "FlowSOM_metacluster", cellular.cols)

##########################################################################################################
#### 8. Annotate clusters
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 5 - annotation")
    setwd("Output 5 - annotation")

    ### Annotate

        annots <- list("T cells" = c(4,1),
                       "Ly6Chi monocyte" = c(12,13,14),
                       "Immature neutrophils" = c(10),
                       "Mature neutriphils" = c(11,15),
                       "B cells" = c(2,3),
                       "Other" = c(5,6,7,9,8)
                       )

        annots <- do.list.switch(annots)
        names(annots) <- c("Values", "Population")
        setorderv(annots, 'Values')
        annots

    ### Add annotations

        cell.dat <- do.add.cols(cell.dat, "FlowSOM_metacluster", annots, "Values")
        cell.dat
        
        cell.sub <- do.add.cols(cell.sub, "FlowSOM_metacluster", annots, "Values")
        cell.sub
        
        make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "Population", col.type = 'factor', add.label = TRUE)
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "Population", group.col, col.type = 'factor')

    ### Expression heatmap
        
        rm(exp)
        exp <- do.aggregate(cell.dat, cellular.cols, by = "Population")
        make.pheatmap(exp, "Population", cellular.cols)
        
    ### Write FCS files

        fwrite(cell.dat, "Annotated.data.csv")
        fwrite(cell.sub, "Annotated.data.DR.csv")
        
        dir.create('FCS files')
        setwd('FCS files')
        
        write.files(cell.dat,
                    file.prefix = exp.name,
                    divide.by = sample.col,
                    write.csv = FALSE,
                    write.fcs = TRUE)

##########################################################################################################
#### Output session info
##########################################################################################################

    ### Session info and metadata
        
        setwd(OutputDirectory)
        dir.create("Output - info", showWarnings = FALSE)
        setwd("Output - info")

        sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
        sessionInfo()
        sink()

        write(raw.cols, "raw.cols.txt")
        write(aligned.cols, "aligned.cols.txt")
        write(cellular.cols, "cellular.cols.txt")
        write(cluster.cols, "cluster.cols.txt")

##########################################################################################################
#### Spectre - Simple Discovery Workflow
##########################################################################################################

    # Spectre R package: https://immunedynamics.io/spectre
    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### Create a folder structure for your analysis run
##########################################################################################################

    ### Create a master folder with a meaningful name. Then inside that folder, insert the following:

        # One folder called 'data' -- this will contain your data CSV or FCS files
        # One folder called 'metadata' -- this will contain a CSV containg your sample metadata
        # One folder called 'Spectre simple discovery' or similar -- place this analysis script there

    ### Example:

        # CNS analysis
        #   /data
        #       -- Contains data files, one CSV or FCS per sample
        #   /metadata
        #       -- Contains a CSV containing sample metadata (group, batch, etc)
        #   /Spectre simple discovery
        #       -- Spectre simple.discovery.R

##########################################################################################################
#### 1. Load packages, and set directories
##########################################################################################################

    ### Load libraries

        library(Spectre)
        Spectre::package.check()    # Check that all required packages are installed
        Spectre::package.load()     # Load required packages

    ### Set PrimaryDirectory
        
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

    ### Set or create 'input' directory
        
        setwd(PrimaryDirectory)
        dir.create('../data', showWarnings = FALSE)
        setwd("../data/")
        InputDirectory <- getwd()
        setwd(PrimaryDirectory)

    ### Set or create 'metadata' directory
        
        setwd(PrimaryDirectory)
        dir.create('../metadata', showWarnings = FALSE)
        setwd("../metadata/")
        MetaDirectory <- getwd()
        setwd(PrimaryDirectory)

    ### Create output directory
        
        setwd(PrimaryDirectory)
        dir.create("Output_Spectre", showWarnings = FALSE)
        setwd("Output_Spectre")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### 2. Import and prep data
##########################################################################################################
    
    ### If you need the demo dataset, uncomment the following code (select all, CMD+SHIFT+C) and run to download.
    ### Alternative: download from https://github.com/ImmuneDynamics/data/blob/main/msCNS.zip?raw=TRUE
        
        # setwd(PrimaryDirectory)
        # setwd("../")
        # getwd()
        # download.file(url = "https://github.com/ImmuneDynamics/data/blob/main/msCNS.zip?raw=TRUE", destfile = 'msCNS.zip', mode = 'wb')
        # unzip(zipfile = 'msCNS.zip')
        # for(i in list.files('msCNS/data', full.names = TRUE)){
        #   file.rename(from = i,  to = gsub('msCNS/', '', i))
        # }
        # for(i in list.files('msCNS/metadata', full.names = TRUE)){
        #   file.rename(from = i,  to = gsub('msCNS/', '', i))
        # }
        # unlink(c('msCNS/', 'msCNS.zip', '__MACOSX'), recursive = TRUE)
        
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

        to.asinh <- names(cell.dat)[c(1:9)]
        to.asinh

        cofactor <- 500
        plot.against <- "Ly6C_asinh"

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
        counts <- meta.dat[,c(2,5)]
        counts

        cell.dat <- do.add.cols(cell.dat, "FileName", sample.info, "Filename", rmv.ext = TRUE)
        cell.dat

    ### Columns

        as.matrix(names(cell.dat))

        cellular.cols <- names(cell.dat)[c(12:20)]
        as.matrix(cellular.cols)

        cluster.cols <- names(cell.dat)[c(12:20)]
        as.matrix(cluster.cols)

        exp.name <- "CNS experiment"
        sample.col <- "Sample"
        group.col <- "Group"
        batch.col <- "Batch"

    ### Subsample targets per group

        data.frame(table(cell.dat[[group.col]])) # Check number of cells per sample.

        unique(cell.dat[[group.col]])
        
        sub.targets <- c(2000, 20000) # target subsample numbers from each group
        sub.targets

##########################################################################################################
#### 5. Clustering and dimensionality reduction
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 2 - clustering")
    setwd("Output 2 - clustering")

    ### Clustering

        cell.dat <- run.flowsom(cell.dat, cluster.cols, meta.k = 8)
        cell.dat
        
    ### Dimensionality reduction

        cell.sub <- do.subsample(cell.dat, sub.targets, group.col)
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
#### 6. Annotate clusters
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 3 - annotation")
    setwd("Output 3 - annotation")

    ### Annotate

        annots <- list("CD4 T cells" = c(2),
                       "CD8 T cells" = c(8,3),
                       "NK cells" = c(1),
                       "Neutrophils" = c(7),
                       "Infil Macrophages" = c(4),
                       "Microglia" = c(5,6)
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
        
        setwd(OutputDirectory)
        setwd("Output 3 - annotation")
        
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
#### 7. Summary data and statistical analysis
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 4 - summary data")
    setwd("Output 4 - summary data")

    ### Setup
    
        variance.test <- 'kruskal.test'
        pairwise.test <- "wilcox.test"
    
        comparisons <- list(c("Mock", "WNV"))
        comparisons
        
        grp.order <- c("Mock", "WNV")
        grp.order
    
    ### Select columns to measure MFI
    
        as.matrix(cellular.cols)
        dyn.cols <- cellular.cols[c(5,8)]
        dyn.cols
    
    ### Create summary tables
    
        sum.dat <- create.sumtable(dat = cell.dat, 
                                   sample.col = sample.col,
                                   pop.col = "Population",
                                   use.cols = dyn.cols, 
                                   annot.cols = c(group.col, batch.col), 
                                   counts = counts)
        
    ### Review summary data
        
        sum.dat
        as.matrix(names(sum.dat))
        
        annot.cols <- c(group.col, batch.col)
        
        plot.cols <- names(sum.dat)[c(4:27)]
        plot.cols
        
    ### Reorder summary data and SAVE
        
        sum.dat <- do.reorder(sum.dat, group.col, grp.order)
        sum.dat[,c(1:3)]
        
        fwrite(sum.dat, 'sum.dat.csv')
        
    ### Autographs

        for(i in plot.cols){
            
            measure <- gsub("\\ --.*", "", i)
            measure
            
            pop <- gsub("^[^--]*.-- ", "", i)
            pop
            
            make.autograph(sum.dat,
                           x.axis = group.col,
                           y.axis = i,
                           y.axis.label = measure,
                           
                           grp.order = grp.order,
                           my_comparisons = comparisons,
                           
                           Variance_test = variance.test,
                           Pairwise_test = pairwise.test,
                           
                           title = pop,
                           subtitle = measure,
                           filename = paste0(i, '.pdf'))
            
        }
        
    ### Create a fold change heatmap
        
        ## Z-score calculation
        sum.dat.z <- do.zscore(sum.dat, plot.cols)
        
        ## Group 
        t.first <- match(grp.order, sum.dat.z[[group.col]])
        t.first <- t.first -1
        t.first
        
        ## Make heatmap
        make.pheatmap(sum.dat.z, 
                      sample.col = sample.col, 
                      plot.cols = paste0(plot.cols, '_zscore'), 
                      is.fold = TRUE, 
                      plot.title = 'Z-score',
                      annot.cols = annot.cols,
                      dendrograms = 'column',
                      row.sep = t.first,
                      cutree_cols = 3)

##########################################################################################################
#### 8. Output session info
##########################################################################################################

    ### Session info and metadata
        
        setwd(OutputDirectory)
        dir.create("Output - info", showWarnings = FALSE)
        setwd("Output - info")

        sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
        sessionInfo()
        sink()

        write(cellular.cols, "cellular.cols.txt")
        write(cluster.cols, "cluster.cols.txt")


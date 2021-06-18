##########################################################################################################
#### Spectre - Multi-level analysis Workflow
##########################################################################################################

    # Spectre R package: https://sydneycytometry.org.au/spectre
    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### 1. Load packages, and set working directory
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

    ### Set 'input' directory
        setwd(PrimaryDirectory)
        setwd("data/")
        InputDirectory <- getwd()
        setwd(PrimaryDirectory)

    ### Set 'metadata' directory
        setwd(PrimaryDirectory)
        setwd("metadata/")
        MetaDirectory <- getwd()
        setwd(PrimaryDirectory)

    ### Create output directory
        dir.create("Output_Spectre", showWarnings = FALSE)
        setwd("Output_Spectre")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

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

        cell.dat <- do.add.cols(cell.dat, "FileName", sample.info, "FileName", rmv.ext = TRUE)
        cell.dat

    ### Columns

        as.matrix(names(cell.dat))

        cellular.cols <- names(cell.dat)[c(11:18)]
        as.matrix(cellular.cols)

        cluster.cols <- names(cell.dat)[c(11:18)]
        as.matrix(cluster.cols)

        exp.name <- "BM experiment"
        sample.col <- "Sample"
        group.col <- "Group"
        batch.col <- "Batch"

    ### Subsample targets per group

        data.frame(table(cell.dat[[group.col]])) # Check number of cells per sample.

        sub.targets <- c(10000, 10000) # target subsample numbers from each group
        sub.targets

##########################################################################################################
#### 5. Clustering and dimensionality reduction
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 2 - clustering")
    setwd("Output 2 - clustering")

    ### Clustering

        cell.dat <- run.flowsom(cell.dat, cluster.cols, meta.k = 20)
        fwrite(cell.dat, "clustered.data.csv")

    ### Dimensionality reduction

        cell.sub <- do.subsample(cell.dat, sub.targets, group.col)
        cell.sub <- run.umap(cell.sub, cluster.cols)

        fwrite(cell.sub, "clustered.data.DR.csv")

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

        annots <- list("CD4 T cells" = c(3),
                       "CD8 T cells" = c(2),
                       "NK cells" = c(1),
                       "Neutrophils" = c(8),
                       "Infil Macrophages" = c(4),
                       "Microglia" = c(5,6,7)
        )

        annots <- do.list.switch(annots)
        names(annots) <- c("Values", "Population")
        setorderv(annots, 'Values')
        annots

    ### Add annotations

        cell.dat <- do.add.cols(cell.dat, "FlowSOM_metacluster", annots, "Values")
        cell.dat
        
        fwrite(cell.dat, "Annotated.data.csv")

        cell.sub <- do.add.cols(cell.sub, "FlowSOM_metacluster", annots, "Values")
        cell.sub
        fwrite(cell.sub, "Annotated.data.DR.csv")

        make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "Population", col.type = 'factor', add.label = TRUE)
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "Population", group.col, col.type = 'factor')

    ### Expression heatmap
        
        rm(exp)
        exp <- do.aggregate(cell.dat, cellular.cols, by = "Population")
        make.pheatmap(exp, "Population", cellular.cols)
        
    ### Write FCS files
        
        setwd(OutputDirectory)
        setwd("Output 3 - annotation")
        
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
        session_info()
        sink()

        write(cellular.cols, "cellular.cols.txt")
        write(cluster.cols, "cluster.cols.txt")


###################################################################################
###  Spectre: spatial 3 - cellular analysis
###################################################################################

    ### Load libraries

        library('Spectre')

        Spectre::package.check(type = 'spatial')
        Spectre::package.load(type = 'spatial')

    ### Set PrimaryDirectory

        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

    ### Set InputDirectory

        setwd(PrimaryDirectory)
        setwd("Output 2 - extract cell data/")
        InputDirectory <- getwd()
        InputDirectory
        
    ### Set metadata directory

        setwd(PrimaryDirectory)
        setwd("metadata")
        MetaDirectory <- getwd()
        MetaDirectory
        
    ### Create output directory

        setwd(PrimaryDirectory)
        dir.create("Output 3 - cellular analysis")
        setwd("Output 3 - cellular analysis")
        OutputDirectory <- getwd()
        OutputDirectory
        
###################################################################################
### Read in spatial data object
###################################################################################

    ### Read in spatial data object

        setwd(InputDirectory)

        list.files(getwd(), '.qs')
        
        spatial.dat <- qread('spatial.dat.qs')
        str(spatial.dat, 3)
        
    ### Metadata
        
        setwd(MetaDirectory)
        
        sample.meta <- fread("ROIs and samples.csv")
        sample.meta

###################################################################################
### Extract cellular data from images using masks and add annotations
###################################################################################

    ### Pull cellular data into a separate data.table

        cell.dat <- do.pull.data(spatial.dat, 'CellData')
        cell.dat

    ### Add sample metadata

        cell.dat <- do.add.cols(cell.dat, base.col = 'ROI', add.dat = sample.meta, add.by = 'ROI')
        cell.dat

    ### Cell type annotations

        types <- list("Cell" = '65534',
                      "Background" = '65535'
                      )

        types <- do.list.switch(types)
        names(types) <- c("Values", "CellType")
        types

        cell.dat <- do.add.cols(cell.dat, base.col = 'cell.type', add.dat = types, add.by = 'Values')
        cell.dat

    ### Region annotations

        regions <- list("Background" = '65533',
                        "White pulp" = '65535',
                        "Red pulp" = "65534")

        regions

        regions <- do.list.switch(regions)
        names(regions) <- c("Values", "Region")
        regions

        cell.dat <- do.add.cols(cell.dat, base.col = 'region', add.dat = regions, add.by = 'Values')
        cell.dat
        
    ### Area calculations
        
        str(spatial.dat, 3)
        
        area.table <- do.calculate.area(spatial.dat, region = 'region')
        area.table
        
        for(i in c(1:length(regions[[1]]))){
            # i <- 1
            nm <- regions[[1]][i]
            trg <- which(names(area.table) == nm)
            names(area.table)[trg] <- regions[[2]][i]
        }
        
        area.table
        
        setwd(OutputDirectory)
        fwrite(area.table, 'area.table.csv')
        
###################################################################################
### Filtering
###################################################################################

    ### Filter (if required - especially for multi-cut)

        as.matrix(unique(cell.dat[['CellType']]))

        cell.dat <- cell.dat[cell.dat[['CellType']] != 'Background',]
        cell.dat

        as.matrix(unique(cell.dat[['CellType']]))

###################################################################################
### Setup columns
###################################################################################

    ### Select 'cellular' columns and perform arcsinh transformation

        as.matrix(names(cell.dat))

        cellular.cols <- names(cell.dat)[c(18:30)]
        as.matrix(cellular.cols)
        
    ### Define key columns
        
        as.matrix(names(cell.dat))
        
        roi.col <- 'ROI'
        sample.col <- 'Sample'
        group.col <- 'Group'
        batch.col <- 'Group'
        
    ### Asinh transformation

        cell.dat <- do.asinh(cell.dat, cellular.cols, cofactor = 1)
        cell.dat

        cellular.cols <- paste0(cellular.cols, "_asinh")
        as.matrix(cellular.cols)

    ### Re-scaling

        cell.dat <- do.rescale(cell.dat, cellular.cols)
        cell.dat

        cellular.cols <- paste0(cellular.cols, "_rescaled")
        as.matrix(cellular.cols)

    ### Define clustering cols

        as.matrix(names(cell.dat))

        cluster.cols <- names(cell.dat)[c(50:58,61:62)]
        as.matrix(cluster.cols)

###################################################################################
### Clustering and dimensionality reduction
###################################################################################

    setwd(OutputDirectory)

    ### Run FlowSOM

        cell.dat <- run.flowsom(cell.dat, cluster.cols, meta.k = 10)
        cell.dat

    ### Run DR

        cell.dat <- run.fitsne(cell.dat, cluster.cols, perplexity = 200)
        cell.dat

    ### Individual plots

        make.colour.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', roi.col, 'factor')
        make.colour.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', sample.col, 'factor')
        make.colour.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', group.col, 'factor')
        make.colour.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', batch.col, 'factor')
        
        make.colour.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', 'CellType', 'factor', add.label = TRUE)
        make.colour.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', 'FlowSOM_metacluster', 'factor', add.label = TRUE)
        
    ### Multi plots

        make.multi.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', plot.by = cellular.cols, figure.title = 'Cellular cols')
        make.multi.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', plot.by = cluster.cols, col.type = 'factor', figure.title = 'Clustering cols')

        make.multi.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', 'FlowSOM_metacluster', 'ROI', col.type = 'factor', figure.title = 'FlowSOM_metacluster by ROI')
        make.multi.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', 'CellType', 'ROI', col.type = 'factor', figure.title = 'CellType by ROI')

    ### Expression heatmap

        exp <- do.aggregate(cell.dat, cellular.cols, 'CellType', 'mean')
        make.pheatmap(exp, 'CellType', cellular.cols)
        rm(exp)

        exp <- do.aggregate(cell.dat, cellular.cols, 'FlowSOM_metacluster', 'mean')
        make.pheatmap(exp, 'FlowSOM_metacluster', cellular.cols)
        rm(exp)

###################################################################################
### Save data
###################################################################################

    setwd(OutputDirectory)
    
    ### qs
    
        qsave(spatial.dat, "spatial.dat.qs")
    
    ### data.table
    
        fwrite(cell.dat, 'cell.dat.csv')
        
    ### FCS files
        
        setwd(OutputDirectory)
        dir.create('FCS')
        setwd('FCS')
        
        write.files(cell.dat, file.prefix = 'IMC_', divide.by = sample.col, write.csv = FALSE, write.fcs = TRUE)
        
###################################################################################
### Spatial plotting of cluster data
###################################################################################

    ### Clusters
    
    setwd(OutputDirectory)
    dir.create('Clusters')
    setwd('Clusters')
    
        for(i in unique(cell.dat$ROI)){
            temp <- cell.dat[cell.dat[['ROI']] == i,]
            make.spatial.plot(spatial.dat, i, 'DNA1_Ir191', mask.outlines = 'cell.mask', temp, cell.col = 'FlowSOM_metacluster', cell.col.type = 'factor')
        }
    
    ### Cell type
    
    setwd(OutputDirectory)
    dir.create('Cell type')
    setwd('Cell type')
    
        for(i in unique(cell.dat$ROI)){
            temp <- cell.dat[cell.dat[['ROI']] == i,]
            make.spatial.plot(spatial.dat, i, 'DNA1_Ir191', mask.outlines = 'cell.mask', temp, cell.col = 'CellType', cell.col.type = 'factor')
        }
        
    ### Folder per 'cluster' -- can compare 'cluster 1' across images
    
    # setwd(OutputDirectory)
    # dir.create('Per cluster')
    # setwd('Per cluster')
    #     
    #     for(a in unique(cell.dat[['FlowSOM_metacluster']])){
    #         
    #         setwd(OutputDirectory)
    #         dir.create('Per cluster')
    #         setwd('Per cluster')
    #         dir.create(paste0('Cluster_', a))
    #         setwd(paste0('Cluster_', a))
    #         
    #         temp.clust <- cell.dat[cell.dat[['FlowSOM_metacluster']] == a,]
    #         
    #         for(i in unique(cell.dat$ROI)){
    #             temp.roi <- temp.clust[temp.clust[['ROI']] == i,]
    #             make.spatial.plot(spatial.dat, i, 'DNA1_Ir191', mask.outlines = 'cell.mask', temp.roi, cell.col = 'FlowSOM_metacluster', cell.col.type = 'factor')
    #         }
    #     }
        
        
    ### Folder per 'cell type' -- can compare 'cluster 1' across images
    
    # setwd(OutputDirectory)
    # dir.create('Per CellType')
    # setwd('Per CellType')
    #     
    #     for(a in unique(cell.dat[['CellType']])){
    #         
    #         setwd(OutputDirectory)
    #         dir.create('Per CellType')
    #         setwd('Per CellType')
    #         dir.create(paste0('CellType_', a))
    #         setwd(paste0('CellType_', a))
    #         
    #         temp.clust <- cell.dat[cell.dat[['CellType']] == a,]
    #         
    #         for(i in unique(cell.dat$ROI)){
    #             temp.roi <- temp.clust[temp.clust[['ROI']] == i,]
    #             make.spatial.plot(spatial.dat, i, 'DNA1_Ir191', mask.outlines = 'cell.mask', temp.roi, cell.col = 'CellType', cell.col.type = 'factor')
    #         }
    #     }
        


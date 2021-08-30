###################################################################################
### Spectre: spatial 2 - cellular analysis
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
        setwd("Output 1 - add masks")
        InputDirectory <- getwd()
        InputDirectory
        
    ### Set metadata directory

        setwd(PrimaryDirectory)
        setwd("metadata")
        MetaDirectory <- getwd()
        MetaDirectory
        
    ### Create output directory

        setwd(PrimaryDirectory)
        dir.create("Output 2 - cellular analysis")
        setwd("Output 2 - cellular analysis")
        OutputDirectory <- getwd()
        OutputDirectory
        
###################################################################################
### Read in spatial data object
###################################################################################

    ### Read in spatial data object

        setwd(InputDirectory)
        setwd("QS file/")
        
        list.files(getwd(), '.qs')
        
        spatial.dat <- qread('spatial.dat.qs')
        str(spatial.dat, 3)
        
        spatial.dat[[1]]$DATA$CellData
        
    ### Metadata
        
        setwd(MetaDirectory)
        
        sample.meta <- fread("ROIs and samples.csv")
        sample.meta

###################################################################################
### Extract cellular data from images using masks and add annotations
###################################################################################

    setwd(OutputDirectory)
        
    ### Pull cellular data into a separate data.table

        cell.dat <- do.pull.data(spatial.dat, 'CellData')
        cell.dat

    ### Add sample metadata

        cell.dat <- do.add.cols(cell.dat, base.col = 'ROI', add.dat = sample.meta, add.by = 'ROI')
        cell.dat
        
    ### Annotate the 'cell.type' masks
        
        make.spatial.plot(spatial.dat, 
                          image.roi = 'ROI002', 
                          image.channel = 'DNA1_Ir191', 
                          cell.dat = 'CellData', 
                          cell.col = 'cell.type', 
                          cell.col.type = 'factor')
        
        
        as.matrix(unique(cell.dat$cell.type))

        cell.type.annot <- list('Cell' = c('65534'),
                                'Non-cell' = c('65535')
                                )
        
        cell.type.annot <- do.list.switch(cell.type.annot)
        names(cell.type.annot) <- c('Values', 'Annotated cell type')
        cell.type.annot
        
        cell.dat <- do.add.cols(cell.dat, 'cell.type', cell.type.annot, 'Values')
        cell.dat
        
    ### Annotate the 'cell.type' masks
        
        make.spatial.plot(spatial.dat, 
                          image.roi = 'ROI002', 
                          image.channel = 'DNA1_Ir191', 
                          cell.dat = 'CellData', 
                          cell.col = 'region', 
                          cell.col.type = 'factor')
        
        as.matrix(unique(cell.dat$region))
        
        region.annot <- list('Background' = c('65533'),
                             'Red pulp' = c('65534'),
                             'White pulp' = c('65535')
                             )
        
        region.annot <- do.list.switch(region.annot)
        names(region.annot) <- c('Values', 'Annotated region')
        region.annot
        
        cell.dat <- do.add.cols(cell.dat, 'region', region.annot, 'Values')
        cell.dat
        
###################################################################################
### Filtering
###################################################################################

    ### Filter
        
        as.matrix(unique(cell.dat$`Annotated cell type`))
        
        cell.dat <- cell.dat[cell.dat[['Annotated cell type']] != 'Non-cell',]
        cell.dat

###################################################################################
### Setup columns
###################################################################################

    ### Select 'cellular' columns and perform arcsinh transformation

        as.matrix(names(cell.dat))

        cellular.cols <- names(cell.dat)[c(8:20)]
        as.matrix(cellular.cols)
        
    ### Define key columns
        
        as.matrix(names(cell.dat))
        
        roi.col <- 'ROI'
        sample.col <- 'Sample'
        group.col <- 'Group'
        batch.col <- 'Batch'
        
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

        cluster.cols <- names(cell.dat)[c(41:47)]
        as.matrix(cluster.cols)

###################################################################################
### Clustering and dimensionality reduction
###################################################################################

    setwd(OutputDirectory)

    ### Run FlowSOM

        cell.dat <- run.flowsom(cell.dat, cluster.cols, meta.k = 15)
        cell.dat

    ### Run DR

        cell.dat <- run.fitsne(cell.dat, cluster.cols, perplexity = 200)
        cell.dat

    ### Individual plots

        make.colour.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', roi.col, 'factor')
        make.colour.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', sample.col, 'factor')
        make.colour.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', group.col, 'factor')
        make.colour.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', batch.col, 'factor')
    
        make.colour.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', 'Annotated cell type', 'factor', add.label = TRUE)
        make.colour.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', 'FlowSOM_metacluster', 'factor', add.label = TRUE)
        
    ### Multi plots

        make.multi.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', plot.by = cellular.cols, figure.title = 'Cellular cols')
        make.multi.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', plot.by = cluster.cols, col.type = 'factor', figure.title = 'Clustering cols')

        make.multi.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', 'Annotated cell type', 'ROI', col.type = 'factor', figure.title = 'FlowSOM_metacluster by ROI')
        make.multi.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', 'FlowSOM_metacluster', 'ROI', col.type = 'factor', figure.title = 'FlowSOM_metacluster by ROI')
        
    ### Expression heatmap

        exp <- do.aggregate(cell.dat, cellular.cols, 'FlowSOM_metacluster', 'mean')
        make.pheatmap(exp, 'FlowSOM_metacluster', cellular.cols)
        rm(exp)

###################################################################################
### Annotate clusters
###################################################################################
        
    ### Cluster annotations
        
        cluster.annots <- list('T cells' = c(12,14,15),
                               'B cells' = c(1,10,5,9))
        
        cluster.annots <- do.list.switch(cluster.annots)
        names(cluster.annots) <- c('Values', 'Annotated metacluster')
        cluster.annots
        
    ### Add annotations, fill in NAs
        
        cell.dat <- do.add.cols(cell.dat, 'FlowSOM_metacluster', cluster.annots, 'Values')
        cell.dat[['Annotated metacluster']][is.na(cell.dat[['Annotated metacluster']])] <- 'Other'
        cell.dat
    
    ### Extra plots
        
        make.colour.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', 'Annotated metacluster', 'factor', add.label = TRUE)
        make.multi.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', 'Annotated metacluster', 'ROI', col.type = 'factor', figure.title = 'Annotated metacluster by ROI')
        
###################################################################################
### Save data
###################################################################################

    setwd(OutputDirectory)
    
    ### data.table
    
        fwrite(cell.dat, 'cell.dat.csv')
        
    ### FCS files
        
        setwd(OutputDirectory)
        dir.create('FCS')
        setwd('FCS')
        
        write.files(cell.dat, file.prefix = 'IMC_All', write.csv = FALSE, write.fcs = TRUE)
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
            make.spatial.plot(spatial.dat, i, 'DNA1_Ir191', mask.outlines = 'cell.mask', temp, cell.col = 'Annotated metacluster', cell.col.type = 'factor')
        }
    
    ### Cell type
    
    setwd(OutputDirectory)
    dir.create('Annotated')
    setwd('Annotated')
    
        for(i in unique(cell.dat$ROI)){
            temp <- cell.dat[cell.dat[['ROI']] == i,]
            make.spatial.plot(spatial.dat, i, 'DNA1_Ir191', mask.outlines = 'cell.mask', temp, cell.col = 'Annotated cell type', cell.col.type = 'factor')
        }

    
    ### Regions
    
    setwd(OutputDirectory)
    dir.create('Regions')
    setwd('Regions')
    
    for(i in unique(cell.dat$ROI)){
        temp <- cell.dat[cell.dat[['ROI']] == i,]
        make.spatial.plot(spatial.dat, i, 'DNA1_Ir191', mask.outlines = 'cell.mask', temp, cell.col = 'Annotated region', cell.col.type = 'factor')
    }
    
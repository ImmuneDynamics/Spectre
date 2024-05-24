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
        setwd("Output 1 - add masks/Data")
        InputDirectory <- getwd()
        InputDirectory
        
    ### Set metadata directory

        setwd(PrimaryDirectory)
        setwd("../metadata")
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
        list.files(getwd(), '.qs')
        
        spatial.dat <- qread('spatial.dat.qs')
        str(spatial.dat, 3)
    
    ### Read in CSV data file
        
        cell.dat <- fread('cell.dat.csv')
        cell.dat    
        
###################################################################################
### Add ROI and sample metadata
###################################################################################

    ### Read in metadata
        
        setwd(MetaDirectory)
        
        sample.meta <- fread("sample.metadata.csv")
        sample.meta
        
    ### Add sample metadata
        
        cell.dat <- do.add.cols(cell.dat, base.col = 'ROI', add.dat = sample.meta, add.by = 'ROI')
        cell.dat

###################################################################################
### Arcsinh transformation
###################################################################################
        
    ### Chose columns to apply arcsinh transformation
        
        as.matrix(names(cell.dat))
        
        to.asinh <- names(cell.dat)[c(6:13)]
        to.asinh
        
    ### Arcsinh
        
        cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = 1)
        cell.dat
        
    ### Re-scaling
        
        cell.dat <- do.rescale(cell.dat, paste0(to.asinh, '_asinh'))
        cell.dat

###################################################################################
### Extract cellular data from images using masks and add annotations
###################################################################################

    ### Define cellular cols
        
        as.matrix(names(cell.dat))
        
        cellular.cols <- names(cell.dat)[c(29:36)]
        cellular.cols
        
    ### Select clustering columns
        
        as.matrix(names(cell.dat))
        
        cluster.cols <- names(cell.dat)[c(29:34)]
        cluster.cols
        
    ### Sample, group, etc columns
        
        as.matrix(names(cell.dat))
        
        roi.col <- 'ROI'
        sample.col <- 'Sample'
        group.col <- 'Group'
        batch.col <- 'Batch'

###################################################################################
### Filtering
###################################################################################

    ### Filter out non-cells or background areas
        
        as.matrix(unique(cell.dat$`Annotated cell type`))
        
        cell.dat <- cell.dat[cell.dat[['Annotated cell type']] != 'Non-cell',]
        cell.dat
        
        as.matrix(unique(cell.dat$`Annotated region`))
        
        cell.dat <- cell.dat[cell.dat[['Annotated region']] != 'Background',]
        cell.dat

    ### Filter by size

        summary(cell.dat$Area)
        
        cell.dat
        cell.dat[cell.dat[['Area']] < 10,]
        
        cell.dat
        cell.dat <- cell.dat[cell.dat[['Area']] > 10,]
        cell.dat
        
    ### Filter plots
        
        setwd(OutputDirectory)
        dir.create('Filter check')
        setwd("Filter check")
        
        as.matrix(names(spatial.dat))
        roi.plot <- names(spatial.dat)[1]
        roi.plot
        
        as.matrix(names(spatial.dat[[1]]@RASTERS))
        
        for(i in cellular.cols){
          make.spatial.plot(spatial.dat, 
                            image.roi =  roi.plot,
                            image.channel = gsub('_asinh_rescaled', '', i), 
                            mask.outlines = 'cell.mask',
                            cell.dat = cell.dat[cell.dat[['ROI']] == roi.plot], 
                            cell.col = 'Annotated cell type', 
                            cell.col.type = 'factor')
        }
        
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
        
    # ### Cluster annotations
    #     
    #     cluster.annots <- list('T cells' = c(12,14,15),
    #                            'B cells' = c(1,10,5,9))
    #     
    #     cluster.annots <- do.list.switch(cluster.annots)
    #     names(cluster.annots) <- c('Values', 'Annotated metacluster')
    #     cluster.annots
    #     
    # ### Add annotations, fill in NAs
    #     
    #     cell.dat <- do.add.cols(cell.dat, 'FlowSOM_metacluster', cluster.annots, 'Values')
    #     cell.dat[['Annotated metacluster']][is.na(cell.dat[['Annotated metacluster']])] <- 'Other'
    #     cell.dat
    # 
    # ### Extra plots
    #     
    #     make.colour.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', 'Annotated metacluster', 'factor', add.label = TRUE)
    #     make.multi.plot(cell.dat, 'FItSNE_X', 'FItSNE_Y', 'Annotated metacluster', 'ROI', col.type = 'factor', figure.title = 'Annotated metacluster by ROI')

###################################################################################
### Save data
###################################################################################

    setwd(OutputDirectory)
    dir.create('Data')
    setwd('Data')
  
    ### Save cellular data
    
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
    
        # setwd(OutputDirectory)
        # dir.create('Annotated clusters')
        # setwd('Annotated clusters')
        # 
        #     for(i in unique(cell.dat$ROI)){
        #         temp <- cell.dat[cell.dat[['ROI']] == i,]
        #         make.spatial.plot(spatial.dat, i, 'DNA1_Ir191', mask.outlines = 'cell.mask', temp, cell.col = 'Annotated metacluster', cell.col.type = 'factor')
        #     }

    ### Cell type
    
        setwd(OutputDirectory)
        dir.create('Annotated cell types')
        setwd('Annotated cell types')
    
            for(i in unique(cell.dat$ROI)){
                temp <- cell.dat[cell.dat[['ROI']] == i,]
                make.spatial.plot(spatial.dat, i, 'DNA1_Ir191', mask.outlines = 'cell.mask', temp, cell.col = 'Annotated cell type', cell.col.type = 'factor')
            }


    ### Regions
    
        setwd(OutputDirectory)
        dir.create('Annotated regions')
        setwd('Annotated regions')
    
        for(i in unique(cell.dat$ROI)){
            temp <- cell.dat[cell.dat[['ROI']] == i,]
            make.spatial.plot(spatial.dat, i, 'DNA1_Ir191', mask.outlines = 'cell.mask', temp, cell.col = 'Annotated region', cell.col.type = 'factor')
        }
    
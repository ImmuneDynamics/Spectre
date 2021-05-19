###################################################################################
### Segmentation 2 - add masks
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
        
    ### Create output directory

        setwd(PrimaryDirectory)
        dir.create("Segmented")
        setwd("Segmented")
        OutputDirectory <- getwd()
        OutputDirectory
        
###################################################################################
### Read in spatial data object
###################################################################################

    ### Read in spatial data object

        setwd(PrimaryDirectory)
        setwd("QS")
        
        list.files(getwd(), '.qs')
        
        spatial.dat <- qread('spatial.dat.qs')
        str(spatial.dat, 3)

###################################################################################
### Read in masks files
###################################################################################

    ### Define cell mask extension for different mask types
        
        setwd(PrimaryDirectory)
        setwd("masks/")
        
        all.masks <- list.files(pattern = '.tif')
        as.matrix(all.masks)
        
    ### Import CELL masks
        
        as.matrix(all.masks)
        cell.mask.ext <- '_Object Identities.tif'
        
        cell.masks <- list.files(pattern = cell.mask.ext)
        as.matrix(cell.masks)
        
        spatial.dat <- do.add.masks(spatial.dat = spatial.dat,
                                    mask.loc = getwd(),
                                    masks = cell.masks,
                                    mask.label = 'cell.mask',
                                    mask.ext = cell.mask.ext)
        
        str(spatial.dat, 3)
        
        spatial.dat[[1]]$RASTERS
        spatial.dat[[1]]$MASKS
        
    ### Import CELL TYPE masks
        
        as.matrix(all.masks)
        cell.type.ext <- '_Object Predictions.tif'
        
        cell.type.masks <- list.files(pattern = cell.type.ext)
        as.matrix(cell.type.masks)

        spatial.dat <- do.add.masks(spatial.dat = spatial.dat,
                                    mask.loc = getwd(),
                                    masks = cell.type.masks,
                                    mask.label = 'cell.type',
                                    mask.ext = cell.type.ext)

        str(spatial.dat, 3)
        
    ### Import REGION masks
        
        as.matrix(all.masks)
        region.ext <- '_Simple Segmentation.tif'
        
        region.masks <- list.files(pattern = region.ext)
        as.matrix(region.masks)
        
        spatial.dat <- do.add.masks(spatial.dat = spatial.dat,
                                    mask.loc = getwd(),
                                    masks = region.masks,
                                    mask.label = 'region',
                                    mask.ext = region.ext)
        
        str(spatial.dat, 3)
        
###################################################################################
### Generate polygons and outlines
###################################################################################

    ### Generate polygons and outlines
        
        spatial.dat <- do.create.outlines(spatial.dat, 'cell.mask')
        spatial.dat <- do.create.outlines(spatial.dat, 'cell.type')
        spatial.dat <- do.create.outlines(spatial.dat, 'region')
        
    ### Checks
        
        str(spatial.dat, 3)
        
        spatial.dat[[1]]$RASTERS
        spatial.dat[[1]]$MASKS
        
        base <- 'DNA1_Ir191'
        mask <- 'cell.mask'
        
    ### Mask plots
        
        setwd(OutputDirectory)
        dir.create('Plots - cell masks')
        setwd('Plots - cell masks')
        
        for(i in names(spatial.dat)){
            make.spatial.plot(spatial.dat = spatial.dat, 
                              image.roi = i, 
                              image.channel = base, 
                              mask.outlines = 'cell.mask')
        }
        
    ### Extract cellular data for each cell mask
        
        spatial.dat <- do.extract(spatial.dat, 'cell.mask')
        str(spatial.dat, 3)
        
    ### Mask plots with data points for cell type
        
        setwd(OutputDirectory)
        dir.create('Plots - cell type masks')
        setwd('Plots - cell type masks')
    
        for(i in names(spatial.dat)){
            make.spatial.plot(spatial.dat = spatial.dat, 
                              image.roi = i, 
                              image.channel = base, 
                              mask.outlines = 'cell.mask', 
                              cell.dat = 'CellData', 
                              cell.col = 'cell.type', 
                              cell.col.type = 'factor')
        }
        
    ### Mask plots with data points for region
        
        setwd(OutputDirectory)
        dir.create('Plots - region masks')
        setwd('Plots - region masks')
        
        for(i in names(spatial.dat)){
            make.spatial.plot(spatial.dat = spatial.dat, 
                              image.roi = i, 
                              image.channel = base, 
                              mask.outlines = 'cell.mask', 
                              cell.dat = 'CellData', 
                              cell.col = 'region', 
                              cell.col.type = 'factor')
        }
       
###################################################################################
### Annotate
###################################################################################       

    ### Annotate cell types
        
        cell.type.annot <- list('Cell' = '257',
                                'Non cell' = '514')
        
        cell.type.annot <- do.list.switch(cell.type.annot)
        names(cell.type.annot) <- c('Values', 'Annotated cell type')
        cell.type.annot
        
        for(i in names(spatial.dat)){
            spatial.dat[[i]]$DATA$CellData <- do.add.cols(spatial.dat[[i]]$DATA$CellData, 'cell.type', cell.type.annot, 'Values')
        }
        
        spatial.dat[[1]]$DATA$CellData
        
    ### Annotate regions
        
        region.annot <- list('White pulp' = '257',
                                'Red pulp' = '514',
                                'Background' = '771')
        
        region.annot <- do.list.switch(region.annot)
        names(region.annot) <- c('Values', 'Annotated region')
        region.annot
        
        for(i in names(spatial.dat)){
            spatial.dat[[i]]$DATA$CellData <- do.add.cols(spatial.dat[[i]]$DATA$CellData, 'region', region.annot, 'Values')
        }
        
        spatial.dat[[1]]$DATA$CellData        
  
###################################################################################
### Area calculations
###################################################################################

    ### Area calculations
        
        str(spatial.dat, 3)
        
        area.table <- do.calculate.area(spatial.dat, region = 'region')
        area.table
        
        for(i in c(1:length(region.annot[[1]]))){
            # i <- 1
            nm <- region.annot[[1]][i]
            trg <- which(names(area.table) == nm)
            names(area.table)[trg] <- region.annot[[2]][i]
        }
        
        area.table
        
        setwd(OutputDirectory)
        fwrite(area.table, 'area.table.csv')
              
###################################################################################
### Save spatial data object as qs file
###################################################################################       
        
    ### Output
        
        setwd(OutputDirectory)
        
        qsave(spatial.dat, "spatial.dat_masks.qs")
        
###################################################################################
### Write FCS files
###################################################################################       

        setwd(OutputDirectory)
        dir.create('FCS files')
        setwd('FCS files')
        
        cell.dat <- do.pull.data(spatial.dat, 'CellData')
        cell.dat
        
        write.files(cell.dat, 'spatial_all', write.csv = FALSE, write.fcs = TRUE)
        write.files(cell.dat, 'spatial', divide.by = 'ROI', write.csv = FALSE, write.fcs = TRUE)
        
        
        
        
        
      
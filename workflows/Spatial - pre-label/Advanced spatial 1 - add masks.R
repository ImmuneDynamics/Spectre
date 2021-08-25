###################################################################################
### Spectre: spatial 1 - add masks and extract cellular data
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
    
    ### Set InputDirectory (ROI TIFFs)
    
    setwd(PrimaryDirectory)
    setwd("data/ROIs/")
    InputDirectory <- getwd()
    InputDirectory
    
    ### Set MaskDirectory (ROI mask TIFFs)
    
    setwd(PrimaryDirectory)
    setwd("data/masks")
    MaskDirectory <- getwd()
    MaskDirectory
    
    ### Create output directory
    
    setwd(PrimaryDirectory)
    dir.create("Output 1 - add masks")
    setwd("Output 1 - add masks")
    OutputDirectory <- getwd()
    OutputDirectory

###################################################################################
### Check ROIs and TIFFs
###################################################################################
    
    ### Initialise the spatial data object with channel TIFF files

    setwd(InputDirectory)

    rois <- list.dirs(full.names = FALSE, recursive = FALSE)
    as.matrix(rois)

    tiff.list <- list()

    ### Check channel names

    for(i in rois){
      setwd(InputDirectory)
      setwd(i)
      tiff.list[[i]] <- list.files(getwd())
    }

    t(as.data.frame(tiff.list))

###################################################################################
### Read in TIFF files and create spatial objects
###################################################################################        
    
    ### Read in ROI channel TIFFs

    setwd(InputDirectory)
    spatial.dat <- read.spatial.files(dir = InputDirectory)

    ### Check results

    str(spatial.dat, 3)
    spatial.dat[[1]]@RASTERS

###################################################################################
### Read in masks files
###################################################################################

    ### Define cell mask extension for different mask types
    
        setwd(MaskDirectory)

        all.masks <- list.files(pattern = '.tif')
        as.matrix(all.masks)
      
        mask.types <- list('cell.mask' = '_ilastik_s2_Object Identities.tif',
                           'cell.type' = '_ilastik_s2_Object Predictions.tif',
                           'regions' = '_ilastik_s2_Simple Segmentation.tif')
        
        mask.types
        
    ### Read in masks
        
        for(i in names(mask.types)){
              spatial.dat <- do.add.masks(dat = spatial.dat, 
                                          mask.dir = MaskDirectory, 
                                          mask.pattern = mask.types[[i]], 
                                          mask.label = i)
        }
        
        str(spatial.dat, 3)
        str(spatial.dat[[1]]@MASKS, 3)
        
###################################################################################
### Rename rasters (if required)
###################################################################################

    ### Check channel names
    
        channel.names <- list()
        
        for(i in names(spatial.dat)){
          channel.names[[i]] <- names(spatial.dat[[i]]@RASTERS)
        }
    
        t(as.data.frame(channel.names))
        
    ### List of corrections (first entry is the 'correct' one)
        
        # corrections <- list(c('CD4','Cd4'),
        #                     c('CD8','CD8a')
        #                     )
        
    ### Replace the 'incorrect' names
        
        # for(i in names(spatial.dat)){
        #   # i <- names(spatial.dat)[[1]]
        #   
        #   for(a in c(1:length(corrections))){
        #     # a <- 1
        #     
        #     trg <- which(names(spatial.dat[[i]]@RASTERS) == corrections[[a]][2])
        #     if(length(trg) != 0){
        #       names(spatial.dat[[i]]@RASTERS)[trg] <- corrections[[a]][1]
        #     }
        #   }
        # }
        
    ### Check channel names
        
        # channel.names <- list()
        # 
        # for(i in names(spatial.dat)){
        #   channel.names[[i]] <- names(spatial.dat[[i]]@RASTERS)
        # }
        # 
        # t(as.data.frame(channel.names))      
        
###################################################################################
### Generate polygons and outlines
###################################################################################

    ### Generate polygons and outlines
        
        spatial.dat <- do.create.outlines(dat = spatial.dat, mask.name = 'cell.mask')
        spatial.dat <- do.create.outlines(dat = spatial.dat, mask.name = 'cell.type')
        spatial.dat <- do.create.outlines(dat = spatial.dat, mask.name = 'regions')
        
    ### Checks
        
        str(spatial.dat, 3)
        str(spatial.dat[[1]]@MASKS, 2)
        
###################################################################################
### Mask QC plots
###################################################################################       
    
    ### Mask plot setup
          
        setwd(OutputDirectory)
        dir.create('Plots - cell masks')
        setwd('Plots - cell masks')
        
        as.matrix(names(spatial.dat[[1]]@RASTERS))
        base <- 'DNA1_Ir191'
        base
        
        as.matrix(names(spatial.dat[[1]]@MASKS))
        mask <- 'cell.mask'
        mask    
        
    ### Create plots
        
        for(i in names(spatial.dat)){
          make.spatial.plot(dat = spatial.dat, 
                            image.roi = i, 
                            image.channel = base, 
                            mask.outlines = mask)
        }
        
###################################################################################
### Calculate cellular data and plot
###################################################################################       
    
    ### Calculate cellular data for each cell mask (this step may take some time)
        
        spatial.dat <- do.extract(spatial.dat, 'cell.mask')
        str(spatial.dat, 3)
        
        spatial.dat[[1]]@DATA
        
    ### Factor plot setup
      
        setwd(OutputDirectory)
        dir.create('Plots - factors')
        setwd('Plots - factors')
        
        plot.rois <- names(spatial.dat)[c(1:2)]
        plot.rois
        
        as.matrix(names(spatial.dat[[1]]@RASTERS))
        base <- 'DNA1_Ir191'
        base
        
        plot.factors <- names(spatial.dat[[1]]@MASKS)[-which('cell.mask' == names(spatial.dat[[1]]@MASKS))]
        plot.factors
        
        plot.exp <- names(spatial.dat[[1]]@RASTERS)
        plot.exp
        
    ### Make factor plots
        
        for(i in plot.rois){
          
          setwd(OutputDirectory)
          setwd('Plots - factors')
          dir.create(i)
          setwd(i)
          
          for(a in plot.factors){
            make.spatial.plot(dat = spatial.dat, 
                              image.roi = i, 
                              image.channel = base, 
                              mask.outlines = mask, 
                              cell.dat = 'CellData', 
                              cell.col = a, 
                              cell.col.type = 'factor',
                              title = paste0(a, ' - ', i))
          }
        }
        
    ### Make exp plots
        
        for(i in plot.rois){
          
          setwd(OutputDirectory)
          setwd('Plots - factors')
          dir.create(i)
          setwd(i)
          
          for(a in plot.exp){
            make.spatial.plot(dat = spatial.dat, 
                              image.roi = i, 
                              image.channel = base, 
                              mask.outlines = mask, 
                              cell.dat = 'CellData', 
                              cell.col = a, 
                              title = paste0(a, ' - ', i))
          }
        }
        
        
###################################################################################
### Extract cellular data and annotate
###################################################################################       

    ### Extract cellular data
        
        cell.dat <- do.pull.data(spatial.dat, 'CellData')
        cell.dat
        
    ### Annotate cell types
        
        cell.type.annot <- list('Cell' = c(65534),
                                'Non-cell' = c(65535)
                                )
        
        cell.type.annot <- do.list.switch(cell.type.annot)
        names(cell.type.annot) <- c('Values', 'Annotated cell type')
        
        cell.dat <- do.add.cols(cell.dat, 'cell.type', cell.type.annot, 'Values')
        cell.dat
       
    ### Annotate regions
        
        region.annot <- list('Background' = c(65533),
                                'Red pulp' = c(65534),
                                'White pup' = c(65535)
                                )
        
        region.annot <- do.list.switch(region.annot)
        names(region.annot) <- c('Values', 'Annotated region')
        
        cell.dat <- do.add.cols(cell.dat, 'regions', region.annot, 'Values')
        cell.dat

###################################################################################
### Area calculation
###################################################################################       

        area.totals <- do.calculate.area(spatial.dat)
        area.totals
        
        area.table <- do.calculate.area(spatial.dat, region = 'regions')
        area.table
        
        for(i in region.annot[['Values']]){
          # i <- 65533
          
          nm <- region.annot[region.annot[['Values']] == i,'Values'][[1]]
          new.nm <- region.annot[region.annot[['Values']] == i,'Annotated region'][[1]]
          
          names(area.table)[which(names(area.table) == i)] <- new.nm
          
        }
        
        area.table
        
###################################################################################
### Save data
###################################################################################       

    ### Output QS and CSV file
        
        setwd(OutputDirectory)
        dir.create('Data')
        setwd('Data')
        
        qsave(spatial.dat, "spatial.dat.qs")
        fwrite(cell.dat, 'cell.dat.csv')
        fwrite(area.totals, 'area.totals.csv')
        fwrite(area.table, 'area.table.csv')

    ### Pull cellular data and write FCS file from each ROI independently
    
        setwd(OutputDirectory)
        dir.create('FCS files')
        setwd('FCS files')
        
        for(i in names(spatial.dat)){
    
          ## Extract data and setup cols
          
              tmp <- list()
              tmp[[i]] <- spatial.dat[[i]]
 
              cell.dat <- do.pull.data(tmp, 'CellData')
              cell.dat <- do.asinh(cell.dat, names(spatial.dat[[i]]@RASTERS), cofactor = 1)
          
          ### Invert y axis
              
              all.neg <- function(test) -1*abs(test)
              
              y_invert <- cell.dat[['y']]
              y_invert <- all.neg(y_invert)
              cell.dat[['y_invert']] <- y_invert
    
          ### Write FCS files  
              
              write.files(cell.dat, i, write.csv = FALSE, write.fcs = TRUE)
              rm(cell.dat)
              rm(i)
        }

    
##########################################################################################################
#### Create a folder structure for your analysis run
##########################################################################################################

    ### Create a master folder with a meaningful name. Then inside that folder, insert the following:
    
        # One folder called 'data' -- this will contain your 'ROIs' folder and 'masks' folder
        # One folder called 'Spectre FlowJo' or similar -- place analysis scripts there
    
    ### Example:
    
        # Spleen analysis
        #   /data
        #       /ROIs -- contains one folder per ROI, filled with TIFF files (one per channel)
        #       /masks -- contains mask files (TIFFs). Can have multiple mask types (e.g. cell mask, region mask)
        #   /Spectre FlowJo
        #       -- TIFF to FCS.R

###################################################################################
### Spectre: TIFF to FCS
###################################################################################

    ### Load libraries
        
        library('Spectre')
        
        Spectre::package.check(type = 'spatial')
        Spectre::package.load(type = 'spatial')
        
    ### Extra packages
        
        if(!require('raster')) {install.packages('raster')}
        if(!require('tiff')) {install.packages('tiff')}
        if(!require('rhdf5')) {BiocManager::install("rhdf5")}
        # if(!require('HDF5Array')) {BiocManager::install("HDF5Array")}
        if(!require('s2')) {install.packages('s2')}
        if(!require('sf')) {install.packages('sf')}
        if(!require('stars')) {install.packages('stars')}
        if(!require('sp')) {install.packages('sp')}
        if(!require('exactextractr')) {install.packages('exactextractr')}
        if(!require('qs')) {install.packages('qs')}
        
        library('raster')
        library('tiff')
        library('rhdf5')
        # library('HDF5Array')
        library('s2')
        library('sf')
        library('stars')
        library('sp')
        library('exactextractr')
        library('qs')
        
    ### Set PrimaryDirectory
        
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory
        
    ### Set InputDirectory (ROI TIFFs)
        
        setwd(PrimaryDirectory)
        dir.create('../data/', showWarnings = FALSE)
        dir.create('../data/ROIs/', showWarnings = FALSE)
        setwd("../data/ROIs/")
        InputDirectory <- getwd()
        InputDirectory
        
    ### Set MaskDirectory (ROI mask TIFFs)
        
        setwd(PrimaryDirectory)
        dir.create('../data/', showWarnings = FALSE)
        dir.create('../data/masks/', showWarnings = FALSE)
        setwd("../data/masks")
        MaskDirectory <- getwd()
        MaskDirectory
        
    ### Create output directory
        
        setwd(PrimaryDirectory)
        dir.create("Output - TIFF to FCS")
        setwd("Output - TIFF to FCS")
        OutputDirectory <- getwd()
        OutputDirectory

###################################################################################
### Check ROIs and TIFFs
###################################################################################
        
    ### If you need the demo dataset, uncomment the following code (select all, CMD+SHIFT+C) and run to download
        
        # setwd(PrimaryDirectory)
        # setwd("../")
        # getwd()
        # download.file(url = "https://github.com/ImmuneDynamics/data/blob/main/spatialFlowjo.zip?raw=TRUE", destfile = 'spatialFlowjo.zip', mode = 'wb')
        # unzip(zipfile = 'spatialFlowjo.zip')
        # for(i in list.files('spatialFlowjo/data', full.names = TRUE)){
        #   file.rename(from = i,  to = gsub('spatialFlowjo/', '', i))
        # }
        # unlink(c('spatialFlowjo/', 'spatialFlowjo.zip', '__MACOSX'), recursive = TRUE)
        
    ### Initialise the spatial data object with channel TIFF files
        
        setwd(InputDirectory)
        
        rois <- list.dirs(full.names = FALSE, recursive = FALSE)
        as.matrix(rois)

    ### Check channel names
        
        tiff.list <- list()
        
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
                           'region' = '_ilastik_s2_Simple Segmentation.tif')
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
        
        for(i in names(spatial.dat[[1]]@MASKS)){
            spatial.dat <- do.create.outlines(dat = spatial.dat, mask.name = i)
        }
        
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
        
        all.dat <- do.pull.data(spatial.dat, 'CellData')
        all.dat

###################################################################################
### Save data
###################################################################################       
        
    ### Output QS and CSV file
        
        setwd(OutputDirectory)
        dir.create('Data')
        setwd('Data')
        
        qs::qsave(spatial.dat, "spatial.dat.qs")
        fwrite(all.dat, 'all.dat.csv')
        
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
                
            ## Invert y axis
            
                all.neg <- function(test) -1*abs(test)
                
                y_invert <- cell.dat[['y']]
                y_invert <- all.neg(y_invert)
                cell.dat[['y_invert']] <- y_invert
            
            ## Write FCS files  
            
                write.files(cell.dat, i, write.csv = FALSE, write.fcs = TRUE)
                rm(cell.dat)
                rm(i)
        }


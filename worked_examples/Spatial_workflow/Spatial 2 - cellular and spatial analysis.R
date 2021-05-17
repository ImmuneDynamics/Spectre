###################################################################################
### Spectre: spatial 2 - add masks and extract cellular data
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
        setwd("Output 1 - Setup/")
        InputDirectory <- getwd()
        InputDirectory
        
    ### Set mask directory
        
        setwd(PrimaryDirectory)
        setwd("masks")
        MaskDirectory <- getwd()
        MaskDirectory
        
    ### Create output directory

        setwd(PrimaryDirectory)
        dir.create("Output 2 - extract cell data")
        setwd("Output 2 - extract cell data")
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
### Read in masks files
###################################################################################

    ### Define cell mask extension for different mask types
        
        setwd(MaskDirectory)
        
        all.masks <- list.files(pattern = '.tif')
        all.masks
        
    ### Import CELL masks
        
        as.matrix(all.masks)
        cell.mask.ext <- '_ilastik_s2_Object Identities.tif'
        
        cell.masks <- list.files(pattern = cell.mask.ext)
        as.matrix(cell.masks)
        
        spatial.dat <- do.add.masks(spatial.dat = spatial.dat,
                                    mask.loc = getwd(),
                                    masks = cell.masks,
                                    mask.label = 'cell.mask',
                                    mask.ext = cell.mask.ext)
        
        str(spatial.dat, 3)
        
    ### Import CELL TYPE masks
        
        as.matrix(all.masks)
        cell.type.ext <- '_ilastik_s2_Object Predictions.tif'
        
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
        region.ext <- '_ilastik_s2_Simple Segmentation.tif'
        
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
        
    ### Extract cellular data for each cell mask
        
        spatial.dat <- do.extract(spatial.dat, 'cell.mask')
        
###################################################################################
### Save spatial data object as qs file
###################################################################################       
        
    ### Output
        
        setwd(OutputDirectory)
        
        qsave(spatial.dat, "spatial.dat.qs")
        
        
      
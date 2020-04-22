

### Load packages

    library(Spectre)
    package.check()
    package.load()

    # package.check(group = "spatial")
    # package.load(group = "spatial")

    library(tiff)
    library(rgeos)
    library(raster)

### Set directories

    dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
    getwd()
    start.dr <- getwd()
    start.dr

    setwd("ROIs/")
    roi.dir <- getwd()

    setwd(start.dr)
    setwd("CPoutput/")
    measure.dir <- getwd()

    setwd(start.dr)

### Read in TIFFs

    setwd(roi.dir)
    rois <- list.dirs(roi.dir, full.names = FALSE, recursive = FALSE)
    rois

    setwd(measure.dir)
    masks <- list.files(measure.dir, ".tiff")
    masks

    cell.dat <- fread("cell.csv")
    image.dat <- fread("Image.csv")

    setwd(start.dr)
    setwd("panel")

    panel.dat <- fread("panel.csv")

    spatial.dat <- read.spatial.files(## ROI tiffs
                                      roi.loc = roi.dir,
                                      rois = rois,

                                      ## Mask tiffs (from CP or Ilastik)
                                      mask.loc = measure.dir,
                                      masks = masks,
                                      mask.ext = "_ilastik_s2_Probabilities_mask.tiff",

                                      ## CP outputs
                                      cell.dat = cell.dat,
                                      image.dat = image.dat,
                                      panel.dat = panel.dat)


    names(spatial.dat)

    for(i in spatial.dat$meta.data$roi.names){
        names(spatial.dat$rois[[i]]$per.cell.means)[10] <- "CD45err"
        names(spatial.dat$rois[[i]]$per.cell.means.filtered)[10] <- "CD45err"
    }

    spatial.dat$rois$`20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac`$per.cell.means

###
    spatial.dat$meta.data

    make.spatial.plot(spatial.dat = spatial.dat,
                      roi = "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac",
                      roi.marker = "CD20_Dy161",
                      cell.dat = "per.cell.means",
                      cell.colour = "CD20")


    cell.dat <- do.extract.cell.dat(spatial.dat, "per.cell.means")

    cell.dat




    make.spatial.plot(spatial.dat = spatial.dat,
                      roi = "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac",
                      roi.marker = "CD20_Dy161",
                      cell.dat = cell.dat[cell.dat[["ImageName"]] == "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac",],
                      cell.colour = "CD20")

    test <- do.label(dat = cell.dat,
                     label.name = "CD20pos",
                     marker.name = "CD20",
                     ">=",
                     20)
    test




    make.spatial.plot(spatial.dat = spatial.dat,
                      roi = "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac",
                      roi.marker = "CD20_Dy161",
                      cell.dat = test[test[["ImageName"]] == "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac",],
                      cell.colour = "CD20pos")
    # throws error, not numeric


    make.spatial.plot(spatial.dat = spatial.dat,
                      roi = "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac",
                      roi.marker = "CD20_Dy161",
                      cell.dat = test[test[["ImageName"]] == "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac",],
                      cell.colour = "CD20pos",
                      cell.colour.type = "factor")






















###

    setwd(start.dr)
    save(spatial.dat, file = "spatial.dat.RData")


###################################################################################################
###################################################################################################
###################################################################################################

    load("spatial.dat.RData")

    spatial.dat$rois
    spatial.dat$rois$`20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac`$
    spatial.dat$rois

    spatial.dat$meta.data
    spatial.dat$per.cell.data
    spatial.dat$cell.dat.means
    spatial.dat$cell.dat.means.filtered


    spatial.dat$cell.dat.means



    spatial.dat$meta.data

    make.spatial.plot(spatial.dat = spatial.dat, # spatial data object

                      roi = "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac", # chosen ROI
                      roi.marker = "CD20_Dy161", # chosen marker for colouration

                      cell.dat = spatial.dat$cell.dat.means[spatial.dat$cell.dat.means[["ImageName"]] == "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac_ilastik_s2_Probabilities_mask.tiff",],
                      cell.x = "X",
                      cell.y = "Y",
                      cell.colour = "CD20")


    temp <- spatial.dat$cell.dat.means[spatial.dat$cell.dat.means[["ImageName"]] == "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac_ilastik_s2_Probabilities_mask.tiff",]

    as.matrix(names(temp))

    temp <- run.flowsom(temp, names(temp)[c(8,9,11,12,13,15)])
    temp <- run.umap(temp, names(temp)[c(8,9,11,12,13,15)])


    make.spatial.plot(spatial.dat = spatial.dat, # spatial data objec
                      roi = "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac", # chosen ROI
                      roi.marker = "CD20_Dy161", # chosen marker for colouration

                      cell.dat = temp,
                      cell.x = "X",
                      cell.y = "Y",
                      cell.colour = "FlowSOM_metacluster")

    make.factor.plot(temp, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", add.label = TRUE)
    make.multi.marker.plot(temp, "UMAP_X", "UMAP_Y", names(temp)[(c(6:18))])







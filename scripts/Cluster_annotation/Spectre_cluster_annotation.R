##########################################################################################################
#### Spectre -- Cluster Annotation
##########################################################################################################

    # Spectre R package: https://sydneycytometry.org.au/spectre
    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### SETUP
##########################################################################################################

    ### Load packages from library
        library(Spectre)
        Spectre::package.check() # --> change so that message at the end is "All required packages have been successfully installed"
        Spectre::package.load() # --> change so that message at the end is "All required packages have been successfully loaded"

        session_info()

    ### Set working directory
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

    ### Determine input directory
        InputDirectory <- PrimaryDirectory
        InputDirectory

    ### Set metadata directory
        setwd(PrimaryDirectory)
        setwd("metadata/")
        MetadataDirectory <- getwd()

    ### Create output directory
        setwd(PrimaryDirectory)
        dir.create("Output_Annotated")
        setwd("Output_Annotated")
        OutputDirectory <- getwd()
        OutputDirectory

        setwd(PrimaryDirectory)

##########################################################################################################
#### LOAD DATA
##########################################################################################################

    ### Read data into R
        setwd(InputDirectory)
        list.files(path = getwd(), pattern = ".csv")

        cell.dat <- fread("Clustered_TAdemo.csv")
        cell.dat

        cell.dat.sub <-fread("DimRed_TAdemo.csv")
        cell.dat.sub

    ### Read in metedata
        setwd(MetadataDirectory)
        list.files(getwd(), ".csv")

        annotations <- fread("annotation.csv")
        annotations

    ### Define existing column name for clusters

        as.matrix(names(cell.dat))

        sample.col <- "Sample"
        group.col <- "Group"

        to.measure <- names(cell.dat)[c(2,4:6,8:9,11:13,16:19,21:30,32)]
        to.annotate <- c(sample.col, group.col)

        as.matrix(names(cell.dat.sub))

        Xdim.name <- "UMAP_X"
        Ydim.name <- "UMAP_Y"

        cluster.name <- "FlowSOM_metacluster"

    ### Create a name for NEW column for populations

        population.name <- "Population"

##########################################################################################################
#### ADD ANNOTATIONS
##########################################################################################################

    ### Embed population name columns in cell.dat
        cell.dat <- Spectre::do.embed.columns(x = cell.dat,
                                              type = "data.table",
                                              base.name = cluster.name,
                                              col.name = population.name,
                                              match.to = as.vector(annotations[,1]),
                                              new.cols = as.vector(annotations[,2]))

        cell.dat

    ### Embed population name columns in cell.dat.sub

        cell.dat.sub <- Spectre::do.embed.columns(x = cell.dat.sub,
                                                  type = "data.table",
                                                  base.name = cluster.name,
                                                  col.name = population.name,
                                                  match.to = as.vector(annotations[,1]),
                                                  new.cols = as.vector(annotations[,2]))

        cell.dat.sub

##########################################################################################################
#### GENERATE SOME NEW DATA (based on population names)
##########################################################################################################

    ### Directories

        setwd(OutputDirectory)
        dir.create("Output-annotated-data", showWarnings = FALSE)
        setwd("Output-annotated-data")

    ### Write 'large' dataset

        Spectre::write.files(dat = cell.dat,
                             file.prefix= paste0("Clustered_annotated"), # required
                             write.csv = TRUE,
                             write.fcs = TRUE)

        Spectre::write.files(dat = cell.dat,
                             file.prefix= paste0("Clustered_annotated"), # required
                             divide.by = sample.col,
                             write.csv = FALSE,
                             write.fcs = TRUE)

    ### Write 'subsample' dataset

        Spectre::write.files(dat = cell.dat.sub,
                             file.prefix = paste0("DimRed_annotated"), # required
                             write.csv = TRUE,
                             write.fcs = TRUE)

        Spectre::write.files(dat = cell.dat.sub,
                             file.prefix = paste0("DimRed_annotated"), # required
                             divide.by = sample.col,
                             write.csv = FALSE,
                             write.fcs = TRUE)

##########################################################################################################
#### Some extra plots
##########################################################################################################

    ### Make plots

        setwd(OutputDirectory)
        dir.create("Output-annotated-plots", showWarnings = FALSE)
        setwd("Output-annotated-plots")

        as.matrix(names(cell.dat.sub))

    ### All data
        make.colour.plot(dat = cell.dat.sub,
                         x.axis = Xdim.name,
                         y.axis = Ydim.name,
                         col.axis = population.name,
                         col.type = 'factor',
                         add.label = TRUE)

    ### Group multi plots
        make.multi.plot(dat = cell.dat.sub,
                        x.axis = Xdim.name,
                        y.axis = Ydim.name,
                        plot.by = population.name,
                        divide.by = group.col,
                        col.type = "factor")

##########################################################################################################
#### Make an expression heatmap
##########################################################################################################

    ### Summary data (cluster x marker)
        Spectre::write.sumtables(dat = cell.dat,
                                 sample.col = sample.col,
                                 pop.col = population.name,
                                 measure.col = to.measure,
                                 annot.col = to.annotate,
                                 group.col = group.col,

                                 do.proportions = FALSE,
                                 cell.counts = NULL,
                                 do.mfi.per.sample = TRUE, ###
                                 do.mfi.per.marker = FALSE)


    ### Make an expression pheatmap
        setwd(OutputDirectory)
        setwd("Output-annotated-plots/SumTable-MFI-PerSample/")

        list.files(getwd(), ".csv")

        to.pheatmap <- read.csv("SumTable-MFI-AllSamples.csv")
        to.pheatmap

    ### Save pheatmap

        setwd(OutputDirectory)
        setwd("Output-annotated-plots")

        make.pheatmap(dat = to.pheatmap,
                      file.name = "Expression by population.png",
                      plot.title = "Expression by population",
                      plot.cols = to.measure,
                      sample.col = population.name)


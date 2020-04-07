##########################################################################################################
#### Spectre -- General Discovery Workflow
#### Part 2/3 - Annotation and re-summarise data
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
        InputDirectory <- paste0(PrimaryDirectory, "/Output_Spectre/Output-data/")
        InputDirectory

    ### Set metadata directory

        setwd(PrimaryDirectory)
        setwd("metadata/")
        MetaDirectory <- getwd()
        MetaDirectory

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

    ### Define some parameters

        as.matrix(names(cell.dat.sub))

        sample.name <- "Sample"
        group.name <- "Group"

        cluster.name <- "FlowSOM_metacluster"

        Xdim.name <- "UMAP_X"
        Ydim.name <- "UMAP_Y"

        to.measure <- names(cell.dat)[c(2,4:6,8:9,11:13,16:19,21:30,32)]
        to.annotate <- names(cell.dat)[c(36:37)]

    ### Read in annotation list

        setwd(MetaDirectory)
        list.files(getwd(), ".csv")

        annotations <- read.csv("annotation.csv")
        annotations

##########################################################################################################
#### ADD ANNOTATIONS
##########################################################################################################


    ### Define name of column for annotated populations

        population.name <- "Population"

    ### Embed population name columns in cell.dat
        cell.dat <- Spectre::do.embed.columns(x = cell.dat,
                                              type = "data.table",
                                              base.name = cluster.name,
                                              col.name = population.name,
                                              match.to = as.numeric(annotations[,1]),
                                              new.cols = as.vector(annotations[,2]))

        cell.dat

    ### Embed population name columns in cell.dat.sub

        cell.dat.sub <- Spectre::do.embed.columns(x = cell.dat.sub,
                                                  type = "data.table",
                                                  base.name = cluster.name,
                                                  col.name = population.name,
                                                  match.to = as.numeric(annotations[,1]),
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
                             divide.by = "Sample",
                             write.csv = FALSE,
                             write.fcs = TRUE)

    ### Write 'subsample' dataset

        Spectre::write.files(dat = cell.dat.sub,
                             file.prefix = paste0("DimRed_annotated"), # required
                             write.csv = TRUE,
                             write.fcs = TRUE)

        Spectre::write.files(dat = cell.dat.sub,
                             file.prefix = paste0("DimRed_annotated"), # required
                             divide.by = "Sample",
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
        make.factor.plot(dat = cell.dat.sub,
                         x.axis = Xdim.name,
                         y.axis = Ydim.name,
                         col.axis = population.name,
                         add.label = TRUE)

        ## Group multi plots
        make.multi.plot(dat = cell.dat.sub,
                        x.axis = Xdim.name,
                        y.axis = Ydim.name,
                        col.axis = population.name,
                        type = "factor",
                        plot.by = group.name,
                        align.xy.by = cell.dat.sub,
                        align.col.by = cell.dat.sub)


    ### Summarise data and make an expression heatmap

        as.matrix(names(cell.dat.sub))

        ## Summary data (cluster x marker)
        Spectre::write.sumtables(x = cell.dat,
                                 sample.col = sample.name,
                                 pop.col = population.name,
                                 measure.col = to.measure,
                                 annot.col = to.annotate,
                                 group.col = group.name,

                                 do.frequencies = FALSE,
                                 cell.counts = NULL,
                                 do.mfi.per.sample = TRUE, ###
                                 do.mfi.per.marker = FALSE)


        ## Make an expression pheatmap
        setwd(OutputDirectory)
        setwd("Output-annotated-plots/SumTable-MFI-PerSample/")
        to.pheatmap <- read.csv("SumTable-MFI-AllSamples.csv")
        to.pheatmap

        make.pheatmap(dat = to.pheatmap,
                      file.name = "Expression by population.png",
                      plot.title = "Expression by population",
                      sample.col = "Population",
                      annot.cols = c(),
                      rmv.cols = c(1))


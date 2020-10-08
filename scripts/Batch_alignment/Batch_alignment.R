##########################################################################################################
#### Batch alignment using Spectre-CytoNorm
##########################################################################################################

##########################################################################################################
#### Setup
##########################################################################################################

    ### Load packages
        library(Spectre)
        Spectre::package.check()
        Spectre::package.load()

        library(CytoNorm)

    ### Set PrimaryDirectory

        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

    ### Set input directory

        setwd("data")
        getwd()
        InputDirectory <- getwd()
        setwd(PrimaryDirectory)

    ### Set metadata directory

        setwd("metadata")
        getwd()
        MetadataDirectory <- getwd()
        setwd(PrimaryDirectory)

    ### Create output directory

        dir.create("output-align")
        setwd("output-align")
        getwd()
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### Import data
##########################################################################################################

    ### Import data

        setwd(InputDirectory)

        list.files(getwd(), ".csv")
        dat.list <- read.files(getwd(), do.embed.file.names = TRUE)
        dat.list

    ### Import metadata

        setwd(MetadataDirectory)

        list.files(getwd(), ".csv")
        meta.dat <- read.csv("sample.details.csv")
        meta.dat




    ### Embed metadata and merge files

        cell.dat <- do.merge.files(dat.list)
        cell.dat

        cell.dat <- do.add.cols(cell.dat, "FileName", meta.dat, "Filename", rmv.ext = TRUE)
        cell.dat

##########################################################################################################
#### Setup columns to align, and specify reference and target data
##########################################################################################################

    ### Specify column denoting samples

        as.matrix(names(cell.dat))

        sample.col <- "Sample"
        batch.col <- "Batch"
        reference.col <- "Reference"

    ### Choose columns to align

        as.matrix(names(cell.dat))

        cluster.cols.nums <- c(1,3:13,15:19,21:22)
        cluster.cols <- names(cell.dat)[cluster.cols.nums] # Channels to align
        cluster.cols

        align.cols.nums <- c(1:22)
        align.cols <- names(cell.dat)[align.cols.nums]
        align.cols

        all.cols.nums <- c(1:22)
        all.cols <- names(cell.dat)[all.cols.nums]
        all.cols

    ### Specify reference data

        ref.dat <- cell.dat[cell.dat[["Reference"]] == TRUE,]
        ref.dat

        setwd(OutputDirectory)
        make.multi.plot(dat = cell.dat,
                        x.axis = "DL800.SA.Bio.Ly6G",
                        y.axis = "BV605.Ly6C",
                        plot.by = batch.col,
                        col.type = "factor",
                        divide.by = sample.col)

    ### Specify target data

        target.dat <- cell.dat
        target.dat

##########################################################################################################
#### Run fsom
##########################################################################################################

    ### Run prep.fsom
        ref.fsom <- do.prep.fsom(dat = ref.dat,
                                 use.cols = cluster.cols,
                                 sample.col = sample.col,
                                 batch.col = batch.col,
                                 xdim = 5,
                                 ydim = 5,
                                 nClus = 15,
                                 scale = FALSE,
                                 seed = 2)

##########################################################################################################
#### Check the clustering to determine it's suitability
##########################################################################################################

    ### Review fsom object version of the output

        ## Main data
        ref.fsom$fsom$FlowSOM

        ## List of reference files and batches
        ref.fsom$fsom$files
        ref.fsom$fsom$batches

        ## Metaclustering
        ref.fsom$fsom$metaclustering

    ### Review data.table  version of the output

        ref.fsom$fsom.dt

    ### Make some plots, to determine suitability of the FlowSOM clustering

        temp <- do.subsample(ref.fsom$fsom.dt, targets = 10000)
        nrow(temp)

        temp <- run.umap(temp, cluster.cols)

        make.colour.plot(temp, "UMAP_X", "UMAP_Y", col.axis = "prep.fsom.metacluster", col.type = 'factor', dot.size = 0.5, title = "Raw - metaclusters", add.label = TRUE)
        make.colour.plot(temp, "UMAP_X", "UMAP_Y", col.axis = "prep.fsom.cluster", col.type = 'factor', dot.size = 0.5, title = "Raw - clusters", add.label = TRUE)
        make.colour.plot(temp, "UMAP_X", "UMAP_Y", col.axis = "FileNo", col.type = 'factor', dot.size = 0.5, title = "Raw - sample")
        make.multi.plot(temp, "UMAP_X", "UMAP_Y",plot.by = all.cols, figure.title = "Raw - markers (all)")


##########################################################################################################
#### Run batch alignment using CytoNorm
##########################################################################################################

    ### Align target data using quantile targets from the fsom file

        new.dat <- do.align(ref.dat = ref.fsom,
                            target.dat = target.dat,

                            ## Columns
                            #sample.col = sample.col, # used to divide 'ref' and 'target' data into 'samples'
                            batch.col = batch.col, # Column denoting batches
                            align.cols = align.cols, # Channels to align

                            ## CytoNorm settings
                            method = "CytoNorm",
                            goal = "mean",
                            nQ = 101,

                            ## Writing result files
                            write.ref.fcs = TRUE,
                            write.target.fcs = TRUE)

##########################################################################################################
#### Check results of alignment in REFERENCE samples
##########################################################################################################

    ### Examine new dataset

        new.dat

    ### Split aligned 'ref' samples and plot to check the quality of alignment

        new.dat.ref <- new.dat[new.dat[[reference.col]] == TRUE,]
        new.dat.ref

        aligned.cluster.cols <- cluster.cols
        for(i in c(1:length(aligned.cluster.cols))){
          a <- aligned.cluster.cols[[i]]
          a <- paste0(a, "_aligned")
          aligned.cluster.cols[[i]] <- a
        }

        aligned.align.cols <- cluster.cols
        for(i in c(1:length(aligned.align.cols))){
          a <- aligned.align.cols[[i]]
          a <- paste0(a, "_aligned")
          aligned.align.cols[[i]] <- a
        }

        aligned.all.cols <- cluster.cols
        for(i in c(1:length(aligned.all.cols))){
          a <- aligned.all.cols[[i]]
          a <- paste0(a, "_aligned")
          aligned.all.cols[[i]] <- a
        }

        aligned.cluster.cols
        aligned.align.cols
        aligned.all.cols

    ### Plots

        new.dat.ref <- do.subsample(new.dat.ref, 10000)
        nrow(new.dat.ref)
        new.dat.ref <- run.umap(new.dat.ref, aligned.align.cols)

        make.colour.plot(new.dat.ref, "UMAP_X", "UMAP_Y", col.axis = sample.col, dot.size = 0.5, title = "Aligned - sample")
        make.multi.plot(new.dat.ref, "UMAP_X", "UMAP_Y", plot.by = aligned.all.cols, figure.title = "Aligned - markers (all)")

        for(i in unique(new.dat.ref$Sample)){
          make.colour.plot(new.dat.ref[new.dat.ref[[sample.col]] == i,], "UMAP_X", "UMAP_Y", col.axis = sample.col, dot.size = 0.5, title = paste0("Aligned_Sample_", i), align.xy.by = new.dat.ref, align.col.by = new.dat.ref)
        }



##########################################################################################################
#### Analyse full dataset
##########################################################################################################

        new.dat <- run.flowsom(new.dat, clustering.cols = aligned.cluster.cols)
        new.dat <- do.subsample(new.dat, 10000)
        new.dat <- run.umap(new.dat, aligned.cluster.cols)

        new.dat <- as.data.table(new.dat) ##########

        ### re-ordered

        make.colour.plot(new.dat, "UMAP_X", "UMAP_Y", col.axis = "FlowSOM_metacluster", title = "All samples - aligned - metaclusters (labelled)", add.label = TRUE)
        make.colour.plot(new.dat, "UMAP_X", "UMAP_Y", col.axis = "FlowSOM_metacluster", title = "All samples - aligned - metaclusters")
        make.multi.plot(dat = new.dat, type = "factor", x.axis = "UMAP_X", y.axis = "UMAP_Y", col.axis = batch.col, plot.by = sample.col, figure.title = "All samples - aligned - samples")



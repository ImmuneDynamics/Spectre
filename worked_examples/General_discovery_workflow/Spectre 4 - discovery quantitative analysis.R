##########################################################################################################
#### Spectre: General Discovery Workflow - (4/4) - Quantitative and Statistical Analysis
##########################################################################################################

    # Spectre R package: https://github.com/ImmuneDynamics/Spectre
    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### Analysis session setup
##########################################################################################################

    ### Load packages

        library(Spectre)
        Spectre::package.check()    # Check that all required packages are installed
        Spectre::package.load()     # Load required packages

    ### Set DT threads

        getDTthreads()

    ### Set primary directory

        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

    ### Set input directory

        setwd(PrimaryDirectory)
        setwd("Output 3 - clustering and DR/Output 3.5 - summary data/")
        InputDirectory <- getwd()
        InputDirectory
        setwd(PrimaryDirectory)

    ### Set metadata directory

        setwd(PrimaryDirectory)
        setwd("metadata/")
        MetaDirectory <- getwd()
        MetaDirectory
        setwd(PrimaryDirectory)

    ### Set output directory

        setwd(PrimaryDirectory)
        dir.create("Output 4 - quantitative analysis", showWarnings = FALSE)
        setwd("Output 4 - quantitative analysis")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### Import data
##########################################################################################################

    setwd(InputDirectory)

    ### Read in files

        list.files(getwd(), ".csv")
        
        sum.dat <- fread('sum.dat.csv')

        as.matrix(names(sum.dat))
        sum.dat    

##########################################################################################################
#### Metadata
##########################################################################################################
        
    setwd(MetaDirectory)
      
    ### Read in files
        
        list.files(getwd(), ".csv")
        
        meta.dat <- fread('sample.details.csv')
        
        as.matrix(names(meta.dat))
        meta.dat          
        
##########################################################################################################
#### Preferences
##########################################################################################################

    ### Define colums

        as.matrix(names(sum.dat))

        sample.col <- "Sample"
        group.col <- "Group"
        batch.col <- "Batch"

        annot.cols <- c(group.col, batch.col)
        annot.cols

        plot.cols <- names(sum.dat)[c(4:28)]
        plot.cols

    ### Comparisons and group order

        variance.test <- 'kruskal.test'
        pairwise.test <- "wilcox.test"

        comparisons <- list(c("Mock", "Virus"))
        comparisons
        
        grp.order <- c("Mock", "Virus")
        grp.order

    ### Sort row order

        sum.dat <- do.reorder(sum.dat, group.col, grp.order)
        
        sum.dat[,c(sample.col, annot.cols),with = FALSE]
        as.matrix(unique(sum.dat[[group.col]]))
        
##########################################################################################################
#### Stats calculations
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 4.1 - summary data tables")
    setwd("Output 4.1 - summary data tables")

    ### Z score
        
        sum.dat.z <- do.zscore(sum.dat, plot.cols)
        sum.dat.z[,plot.cols] <- NULL
        names(sum.dat.z) <- gsub('_zscore', '', names(sum.dat.z))
        
        as.matrix(names(sum.dat.z))
        sum.dat.z

        fwrite(sum.dat, "Summary data.csv")
        fwrite(sum.dat.z, "Summary data - z-score.csv")

    ### Statistical tests

        sum.dat.stats.raw <- create.stats(sum.dat,
                                         use.cols = plot.cols,
                                         sample.col = sample.col,
                                         group.col = group.col,
                                         comparisons = comparisons,
                                         corrections = NULL,
                                         variance.test = variance.test,
                                         pairwise.test = pairwise.test)
        
        sum.dat.stats.raw

        
        sum.dat.stats.FDR <- create.stats(sum.dat,
                                        use.cols = plot.cols,
                                        sample.col = sample.col,
                                        group.col = group.col,
                                        comparisons = comparisons,
                                        corrections = 'fdr',
                                        variance.test = variance.test,
                                        pairwise.test = pairwise.test)
        
        sum.dat.stats.FDR

    ### Review and save to disk
        
        sum.dat.stats.raw[,c(1:3)]
        sum.dat.stats.FDR[,c(1:3)]

        fwrite(sum.dat.stats.raw, "Summary data - stats - uncorrected.csv")
        fwrite(sum.dat.stats.FDR, "Summary data - stats - FDR.csv")
        
    ###  p-value tables

        raw <- sum.dat.stats.raw[sum.dat.stats.raw[["Type"]] == 'p-value',]
        raw <- raw[raw[["Comparison"]] != 'Kruskal',]

        pval.compars <- raw[["Comparison"]]

        pval <- raw[,plot.cols, with = FALSE]
        pval.sig <- matrix(nrow = 0, ncol = length(plot.cols))

        for(i in c(1:nrow(pval))){
          temp <- pval[i,]
          temp <- temp < 0.05
          temp <- gsub(TRUE, "Significant", temp)
          temp <- gsub(FALSE, "NS", temp)
          pval.sig <- rbind(pval.sig, temp)
        }

        pval.sig <- as.data.frame(pval.sig)
        names(pval.sig) <- plot.cols
        rownames(pval.sig) <- paste0("p-value - ", pval.compars)

        p.val.annots <- list()

        for(i in rownames(pval.sig)){
          p.val.annots[[i]]  <- c('NS' = "Black", "Significant" = "Blue")
        }

    ### q-value tables (FDR adjusted p-values)

        fdr <- sum.dat.stats.FDR[sum.dat.stats.FDR[["Type"]] == 'p-value_fdr',]
        fdr <- fdr[fdr[["Comparison"]] != 'Kruskal',]

        pval.fdr.compars <- fdr[["Comparison"]]

        pval.fdr <- fdr[,plot.cols, with = FALSE]
        pval.fdr.sig <- matrix(nrow = 0, ncol = length(plot.cols))

        for(i in c(1:nrow(pval.fdr))){
          temp <- pval.fdr[i,]
          temp <- temp < 0.05
          temp <- gsub(TRUE, "Significant", temp)
          temp <- gsub(FALSE, "NS", temp)
          pval.fdr.sig <- rbind(pval.fdr.sig, temp)
        }

        pval.fdr.sig <- as.data.frame(pval.fdr.sig)
        names(pval.fdr.sig) <- plot.cols
        rownames(pval.fdr.sig) <- paste0("q-value - ", pval.fdr.compars)

        p.val.fdr.annots <- list()

        for(i in rownames(pval.fdr.sig)){
          p.val.fdr.annots[[i]]  <- c('NS' = "Black", "Significant" = "Red")
        }

    ### Create annotation data.frame

        feature.annots <- rbind(pval.sig, pval.fdr.sig)
        feature.annots <- t(feature.annots)
        feature.annots <- as.data.frame(feature.annots)
        
        feature.annots
        str(feature.annots)

        p.val.annots
        p.val.fdr.annots

        annotation_colors <- c(p.val.annots, p.val.fdr.annots)

##########################################################################################################
#### Differential heatmap
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 4.2 - heatmaps")
    setwd("Output 4.2 - heatmaps")

    ### Determine rows for the start of each new group
    
        t.first <- match(grp.order, sum.dat.z[[group.col]])
        t.first <- t.first -1
    
    ### Heatmap (static rows)
    
        make.pheatmap(sum.dat.z,
                      sample.col = sample.col,
                      plot.cols = plot.cols,
                      annot.cols = annot.cols,
                      feature.annots = feature.annots,
                      annotation_colors = annotation_colors,
                      is.fold = TRUE,
                      fold.range = c(3, -3),
                      dendrograms = 'column',
                      row.sep = t.first,
                      #cutree_rows = 2,
                      cutree_cols = 3,
                      plot.title = "All features - z-score (static rows)",
                      file.name = "All features - z-score (static rows).png")
        
    ### Heatmap (clustered rows)
        
        make.pheatmap(sum.dat.z,
                      sample.col = sample.col,
                      plot.cols = plot.cols,
                      annot.cols = annot.cols,
                      feature.annots = feature.annots,
                      annotation_colors = annotation_colors,
                      is.fold = TRUE,
                      fold.range = c(3, -3),
                      cutree_rows = 2,
                      cutree_cols = 3,
                      plot.title = "All features - z-score",
                      file.name = "All features - z-score.png")

##########################################################################################################
#### PCA
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 4.3 - PCA")
    setwd("Output 4.3 - PCA")
    
    ### Setup
    
        any(is.na(sum.dat))
    
    ### PCA plots - columns with NA removed
        
        setwd(OutputDirectory)
        setwd("Output 4.3 - PCA")
        dir.create('PCA -- NA cols removed')
        setwd('PCA -- NA cols removed')
        
        for.pca <- sum.dat[, plot.cols, with = FALSE]
        for.pca <- for.pca[ , colSums(is.na(for.pca)) == 0 , with = FALSE]
        
        pca.cols <- names(for.pca)
        pca.cols
        
        for.pca <- cbind(sum.dat[,c(sample.col, annot.cols), with = FALSE], for.pca)
        for.pca

        pca_out <- stats::prcomp(for.pca[,pca.cols,with = FALSE],
                                 scale = TRUE)
        
        pcs <- colnames(pca_out$x)[c(1:length(colnames(pca_out$x)))]
        pcs
        
        pca.plotting <- cbind(for.pca, pca_out$x[,pcs])
        
        make.colour.plot(pca.plotting, 'PC1', 'PC2', group.col, 'factor', dot.size = 4)
        
        rm(for.pca)
        rm(pca.cols)
        rm(pca_out)
        rm(pcs)
        rm(pca.plotting)
        
        # Spectre::run.pca(dat = for.pca,
        #                  use.cols = plot.cols,
        #                  # scale = FALSE,
        #                  # plot.ind.label = c("point", "text"),
        #                  # row.names = patient.id,
        #                  plot.ind.group = TRUE,
        #                  group.ind = group.col,
        #                  repel = FALSE)
        
    ### PCA plots - rows with NA removed
    
        setwd(OutputDirectory)
        setwd("Output 4.3 - PCA")
        dir.create('PCA -- NA rows removed')
        setwd('PCA -- NA rows removed')
    
        for.pca <- sum.dat[, plot.cols, with = FALSE]
        cmpl.rws <- complete.cases(for.pca)
        for.pca <- for.pca[cmpl.rws, ]
    
        pca.cols <- names(for.pca)
        pca.cols
        
        for.pca <- cbind(sum.dat[,c(sample.col, annot.cols), with = FALSE], for.pca)
        for.pca
    
        pca_out <- stats::prcomp(for.pca[,pca.cols,with = FALSE],
                                 scale = TRUE)
        
        pcs <- colnames(pca_out$x)[c(1:length(colnames(pca_out$x)))]
        pcs
        
        pca.plotting <- cbind(for.pca, pca_out$x[,pcs])
        
        make.colour.plot(pca.plotting, 'PC1', 'PC2', group.col, 'factor', dot.size = 4)
        
        rm(for.pca)
        rm(pca.cols)
        rm(pca_out)
        rm(pcs)
        rm(pca.plotting)
        
        # Spectre::run.pca(dat = for.pca,
        #                  use.cols = plot.cols,
        #                  # scale = FALSE,
        #                  # plot.ind.label = c("point", "text"),
        #                  # row.names = patient.id,
        #                  plot.ind.group = TRUE,
        #                  group.ind = group.col,
        #                  repel = FALSE)

##########################################################################################################
#### Volcano plots
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 4.4 - volcano plots")
    setwd("Output 4.4 - volcano plots")

    ### Setup for volcanos

        comps <- list()

        for(i in c(1:length(comparisons))){
          temp <- comparisons[[i]]
          strg <- paste0(temp[[1]], " to ", temp[[2]])
          comps[[i]] <- strg
        }

        comps

    ### Uncorrected volcanos (p-values)
        
        setwd(OutputDirectory)
        dir.create("Output 4.4 - volcano plots")
        setwd("Output 4.4 - volcano plots")

        dir.create("Uncorrected p-values")
        setwd("Uncorrected p-values")

        for(i in comps){

          temp <- sum.dat.stats.raw[sum.dat.stats.raw[["Comparison"]] == i,]

          p.dat <- temp[temp[["Type"]] == "p-value",]
          p.dat <- p.dat[,names(p.dat)[c(3:length(names(p.dat)))], with = FALSE]

          fc.dat <- temp[temp[["Type"]] == "FClog2",]
          fc.dat <- fc.dat[,names(fc.dat)[c(3:length(names(fc.dat)))], with = FALSE]

          nms <- names(fc.dat)

          make.volcano.plot(dat.p = p.dat,
                            dat.fc = fc.dat,
                            vars = nms,
                            title = i,
                            xlim = c(-3.5, 3.5))
        }

    ### Corrected volcanos (q-values)
        
        setwd(OutputDirectory)
        dir.create("Output 4.4 - volcano plots")
        setwd("Output 4.4 - volcano plots")

        dir.create("Corrected p-values")
        setwd("Corrected p-values")

        for(i in comps){

          temp <- sum.dat.stats.FDR[sum.dat.stats.FDR[["Comparison"]] == i,]

          p.dat <- temp[temp[["Type"]] == "p-value_fdr",]
          p.dat <- p.dat[,names(p.dat)[c(3:length(names(p.dat)))], with = FALSE]

          fc.dat <- temp[temp[["Type"]] == "FClog2",]
          fc.dat <- fc.dat[,names(fc.dat)[c(3:length(names(fc.dat)))], with = FALSE]

          nms <- names(fc.dat)

          make.volcano.plot(dat.p = p.dat,
                            dat.fc = fc.dat,
                            vars = nms,
                            title = i,
                            xlim = c(-3.5, 3.5))
        }

##########################################################################################################
#### Violin plots (with stats)
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 4.5 - autographs")
    setwd("Output 4.5 - autographs")

    # meas.type <- unique(sub(" -- .*", "", names(sum.dat[,..plot.cols])))
    # meas.type

    for(i in plot.cols){

      pop <- sub(".* -- ", "", i) # population
      meas <- sub(" -- .*", "", i) # measurement

      make.autograph(sum.dat,
                     x.axis = group.col,
                     y.axis = i,
                     y.axis.label = meas,

                     grp.order = grp.order,
                     my_comparisons = comparisons,

                     Variance_test = variance.test,
                     Pairwise_test = pairwise.test,

                     title = pop,
                     subtitle = meas
      )
    }

##########################################################################################################
#### Output session info
##########################################################################################################

    ### Save session info
        
        setwd(OutputDirectory)
        dir.create("Output - info")
        setwd("Output - info")
        
        sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
        session_info()
        sink()
        
        setwd(PrimaryDirectory)
    
    
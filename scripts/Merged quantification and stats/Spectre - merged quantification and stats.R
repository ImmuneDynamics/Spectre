##########################################################################################################
#### Spectre discovery workflow - C - differential quantiative and statistical analysis
##########################################################################################################

    # Spectre R package: https://sydneycytometry.org.au/spectre
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
        setwd("Output B - discovery analysis/Output - summary data/")
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
        dir.create("Output C - quantitative analysis", showWarnings = FALSE)
        setwd("Output C - quantitative analysis")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

    ### Save session info

        setwd(OutputDirectory)
        dir.create("Output - info")
        setwd("Output - info")

        sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
        session_info()
        sink()

        setwd(PrimaryDirectory)

##########################################################################################################
#### Import data
##########################################################################################################

    setwd(InputDirectory)

    ### Read in files
        sum.list <- list()

        files <- list.files(getwd(), ".csv")
        files

        for(i in files){
          sum.list[[i]] <- fread(i)
        }

        sum.list[[1]]

    ### Merge into single data.table

        as.matrix(names(sum.list[[1]]))

        checks <-do.list.summary(sum.list)
        checks$name.table

        metric.col <- c(1)
        essential.cols <- c(2,3,4)
        rmv.cols <- c(5)

        sum.dat <- sum.list[[1]][,c(essential.cols), with = FALSE]
        sum.dat

        for(i in files){
          temp <- sum.list[[i]]
          metric <- temp[1,metric.col, with = FALSE]

          temp <- temp[,-c(metric.col, essential.cols, rmv.cols), with = FALSE]
          names(temp) <- paste0(metric, " -- ", names(temp))

          sum.dat <- cbind(sum.dat, temp)
        }

    ### Check reults

        sum.dat
        as.matrix(names(sum.dat))

##########################################################################################################
#### Preferences
##########################################################################################################

    ### Preferences

        as.matrix(names(sum.dat))

        sample.col <- "Sample"
        group.col <- "Group"
        batch.col <- "Batch"

        annot.cols <- c(group.col, batch.col)
        annot.cols

        plot.cols <- names(sum.dat)[c(4:9,64:69)]
        plot.cols

    ### Comparisons

        variance.test <- NULL # 'kruskal.test'
        pairwise.test <- "wilcox.test"

        comparisons <- list(c("Mock", "WNV"))

        grp.order <- c("Mock", "WNV")

    ### Other adjustments

        sum.dat$Batch <- as.factor(sum.dat$Batch)

##########################################################################################################
#### Stats calculations
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output - summary data tables")
    setwd("Output - summary data tables")

    ### Z score

        zscore <- function(x) {
          z <- (x - mean(x)) / sd(x)
          return(z)
        }

        sum.dat.z <- sum.dat

        res <- scale(sum.dat.z[,plot.cols, with = FALSE])
        res <- as.data.table(res)
        sum.dat.z[,plot.cols] <- res

        sum.dat.z

        fwrite(sum.dat, "Summary data.csv")
        fwrite(sum.dat.z, "Summary data - z-score.csv")

    ### Statistical tests

        sum.dat.stats.raw <- do.stats(sum.dat,
                                       use.cols = plot.cols,
                                       sample.col = sample.col,
                                       grp.col = group.col,
                                       comparisons = comparisons,
                                       corrections = NULL,
                                       variance.test = variance.test,
                                       pairwise.test = pairwise.test)

        sum.dat.stats.FDR <- do.stats(sum.dat,
                                      use.cols = plot.cols,
                                      sample.col = sample.col,
                                      grp.col = group.col,
                                      comparisons = comparisons,
                                      corrections = 'fdr',
                                      variance.test = variance.test,
                                      pairwise.test = pairwise.test)
        sum.dat.stats.raw
        sum.dat.stats.FDR

        sum.dat.stats.raw[,c(1:3)]
        sum.dat.stats.FDR[,c(1:3)]

        fwrite(sum.dat.stats.raw, "Summary data - stats -  uncorrected.csv")
        fwrite(sum.dat.stats.FDR, "Summary data - stats - FDR.csv")

    ### Sig tables

        ## Raw p-values

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

        ## P-values FDR

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
            rownames(pval.fdr.sig) <- paste0("p-value_fdr - ", pval.fdr.compars)

            p.val.fdr.annots <- list()

            for(i in rownames(pval.fdr.sig)){
              p.val.fdr.annots[[i]]  <- c('NS' = "Black", "Significant" = "Red")
            }


        ## Create annotation data.frame

            x <- rbind(pval.sig, pval.fdr.sig)

            # x <- data.frame("p_value" = pval,
            #                 "p_value_FDR" = pval.fdr)

            x <- t(x)
            x <- as.data.frame(x)

            # x <- as.matrix(x)
            # feature.annots <- as.data.frame(x)
            # rownames(feature.annots) <- plot.cols
            # feature.annots

            #str(feature.annots)
            #str(my_sample_col)
            str(x)

            feature.annots <- x

            p.val.annots
            p.val.fdr.annots

            annotation_colors <- c(p.val.annots, p.val.fdr.annots)

            # annotation_colors <- list('p_value' = c('NS' = "Black", "Significant" = "Blue"),
            #                           'p_value_FDR' = c('NS' = "Black", "Significant" = "Red"))


##########################################################################################################
#### Differential heatmap
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output - heatmaps")
    setwd("Output - heatmaps")

    sum.dat.z[[group.col]]

    make.pheatmap(sum.dat.z,
                  sample.col = sample.col,
                  plot.cols = plot.cols,
                  annot.cols = annot.cols,
                  feature.annots = feature.annots,
                  annotation_colors = annotation_colors,
                  is.fold = TRUE,
                  fold.range = c(3, -3),
                  dendrograms = 'column',
                  row.sep = 6,
                  #cutree_rows = 2,
                  cutree_cols = 2,
                  plot.title = "All features - z-score (static rows)",
                  file.name = "All features - z-score (static rows).png")

    make.pheatmap(sum.dat.z,
                  sample.col = sample.col,
                  plot.cols = plot.cols,
                  annot.cols = annot.cols,
                  feature.annots = feature.annots,
                  annotation_colors = annotation_colors,
                  is.fold = TRUE,
                  fold.range = c(3, -3),
                  cutree_rows = 2,
                  cutree_cols = 2,
                  plot.title = "All features - z-score",
                  file.name = "All features - z-score.png")

##########################################################################################################
#### PCA
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output - PCA")
    setwd("Output - PCA")

    ### Select for relevant columns
        for.pca <- sum.dat[, c(sample.col, annot.cols, plot.cols), with = FALSE]
        for.pca

        any(is.na(for.pca))

    ### Remove NAs

        ## Remove columns with NA
            # for.pca <- for.pca[ , colSums(is.na(for.pca)) == 0 , with = FALSE]
            # pca.cols <- names(for.pca)[c(3:length(names(for.pca)))]

        ## Remove ROWS with NA
            # for.pca <- for.pca[complete.cases(for.pca), ]
            # pca.cols <- names(for.pca)[c(4:length(names(for.pca)))]

        any(is.na(for.pca))

    ### Simple PCA plot

        pca_out <- stats::prcomp(for.pca[,plot.cols,with = FALSE],
                                 scale = TRUE)

        pcs <- colnames(pca_out$x)[c(1:10)]
        pcs

        pca.plotting <- cbind(for.pca, pca_out$x[,c(1:10)])
        make.colour.plot(pca.plotting, 'PC1', 'PC2', group.col, dot.size = 5)

    ### Complex PCA
        Spectre::run.pca(dat = for.pca,
                         use.cols = plot.cols,
                         # scale = FALSE,
                         # plot.ind.label = c("point", "text"),
                         # row.names = patient.id,
                         plot.ind.group = TRUE,
                         group.ind = group.col,
                         repel = FALSE)

##########################################################################################################
#### Volcano plots
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output - volcano plots")
    setwd("Output - volcano plots")

    ### Setup

        comps <- list()

        for(i in c(1:length(comparisons))){
          temp <- comparisons[[i]]
          strg <- paste0(temp[[1]], " to ", temp[[2]])
          comps[[i]] <- strg
        }

        comps

    ### Uncorrected volcanos
        setwd(OutputDirectory)
        dir.create("Output - volcano plots")
        setwd("Output - volcano plots")

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

    ### Corrected
        setwd(OutputDirectory)
        dir.create("Output - volcano plots")
        setwd("Output - volcano plots")

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
#### AutoGraphs
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output - autographs")
    setwd("Output - autographs")

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

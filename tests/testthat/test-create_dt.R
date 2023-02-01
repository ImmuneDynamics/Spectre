#' Unit test for creating data.table from other popular SC based objects.

# Install from using devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(Seurat)

library(SingleCellExperiment)
library(scRNAseq)
library(Rtsne)

test_that("Non supported objects are not converted", {
    to_convert <- c(1:20)
    expect_error(create.dt(to_convert, from = "Vector"))
    expect_error(create.dt(to_convert))
})

test_that("Seurat RNA only object is converted", {
    # Ok suppress warnings is not great, but seurat data need to fix this!
    suppressWarnings(InstallData("pbmc3k"), classes = "message")
    pbmc <- suppressWarnings(LoadData("pbmc3k", type = "pbmc3k.final"), classes = "message")

    # Should succeed
    expect_error(create.dt(pbmc, from = "Seurat"), NA)

    dat <- suppressWarnings(create.dt(pbmc), classes = "message")

    # Meant to return a list
    expect_equal(typeof(dat), "list")

    # Check all the genes are captured
    gene_names <- paste0(rownames(pbmc), "_RNA_counts")
    expect_true(all(gene_names %in% colnames(dat[["data.table"]])))
    # There is an element in the list that stores the gene names
    expect_true(all(rownames(pbmc) %in% dat[["geneNames"]]))

    # Check all cell barcodes are captured
    expect_true(all(colnames(pbmc) %in% dat[["data.table"]][["cellNames"]]))
    expect_true(all(colnames(pbmc) %in% dat[["cellNames"]]))

    # Check metadata
    expect_true(all(colnames(pbmc@meta.data) %in% dat[["meta.data"]]))
    for (meta_dat in colnames(pbmc@meta.data)) {
        expected <- pbmc@meta.data[[meta_dat]]
        obj <- dat[["data.table"]][[meta_dat]]
        expect_true(all.equal(obj, expected))
    }

    # Check HVG
    expect_true(all(VariableFeatures(pbmc) %in% dat[["var.features"]]))
    expect_true(all(head(VariableFeatures(pbmc), 10) %in% dat[["var.features.top10"]]))

    # Check assay
    expect_true("_RNA" == dat[["assays"]])

    # Check all slots are captured
    expect_true(all(c("_counts", "_data", "_scale.data") %in% dat[["slots"]]))

    # Check dimreds
    expect_true(all(c("pca_", "umap_") %in% dat[["dim.reds"]]))
    # By default the PBMC data has 50 PCs.
    pca_coordinates <- Embeddings(pbmc, reduction = "pca")
    for (i in c(1:50)) {
        expected <- pca_coordinates[, paste0("PC_", i)]
        obj <- dat[["data.table"]][[paste0("pca_PC_", i)]]
        names(obj) <- dat[["data.table"]][["cellNames"]]
        expect_true(identical(expected, obj))
    }

    # By default, the PBMC data has 2 UMAP
    umap_coordinates <- Embeddings(pbmc, reduction = "umap")
    for (i in c(1:2)) {
        expected <- umap_coordinates[, paste0("UMAP_", i)]
        obj <- dat[["data.table"]][[paste0("umap_UMAP_", i)]]
        names(obj) <- dat[["data.table"]][["cellNames"]]
        expect_true(identical(expected, obj))
    }
})

test_that("SCE object is converted", {
    sce <- suppressWarnings(ReprocessedAllenData("tophat_counts"), classes = "message")
    # Quick prep to get counts and dim red.
    # See https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html
    counts(sce) <- assay(sce, "tophat_counts")
    assay(sce, "tophat_counts") <- NULL

    libsizes <- colSums(counts(sce))
    size.factors <- libsizes / mean(libsizes)
    logcounts(sce) <- log2(t(t(counts(sce)) / size.factors) + 1)
    pca_data <- prcomp(t(logcounts(sce)), rank = 50)
    set.seed(42)
    tsne_data <- Rtsne(pca_data$x[, 1:50], pca = FALSE)
    reducedDims(sce) <- list(PCA = pca_data$x, TSNE = tsne_data$Y)

    # Should succeed
    expect_error(create.dt(sce, from = "SingleCellExperiment"), NA)

    dat <- suppressWarnings(create.dt(sce), classes = "message")

    # Meant to return a list
    expect_equal(typeof(dat), "list")

    # Check number of cells
    expect_equal(dim(assay(sce))[2], nrow(dat[["data.table"]]))

    # Check all the genes are captured
    gene_names <- c(
        paste0(rownames(sce), "_counts"),
        paste0(rownames(sce), "_logcounts")
    )
    expect_true(all(gene_names %in% colnames(dat[["data.table"]])))
    # There is an element in the list that stores the gene names
    expect_true(all(rownames(sce) %in% dat[["geneNames"]]))

    # Check all cell barcodes are captured
    expect_true(all(colnames(sce) %in% dat[["data.table"]][["cellNames"]]))
    expect_true(all(colnames(sce) %in% dat[["cellNames"]]))

    # Check metadata
    expect_true(all(colnames(colData(sce)) %in% dat[["meta.data"]]))
    for (meta_dat in colnames(colData(sce))) {
        expected <- colData(sce)[[meta_dat]]
        obj <- dat[["data.table"]][[meta_dat]]
        expect_true(all.equal(obj, expected))
    }

    # Check all counts assays are captured
    expect_true(all(c("_counts", "_logcounts") %in% dat[["assays"]]))

    # Check dimreds
    expect_true(all(c("PCA_", "TSNE_") %in% dat[["dim.reds"]]))

    # We set up the PCA to have 50 PCs
    pca_coordinates <- reducedDim(sce, "PCA")
    for (i in c(1:50)) {
        expected <- pca_coordinates[, paste0("PC", i)]
        obj <- dat[["data.table"]][[paste0("PCA_", i)]]
        names(obj) <- dat[["data.table"]][["cellNames"]]
        expect_true(identical(expected, obj))
    }

    # By default, the data has 2 TSNE coordinates
    tsne_coordinates <- reducedDim(sce, "TSNE")
    for (i in c(1:2)) {
        expected <- tsne_coordinates[, i]
        obj <- dat[["data.table"]][[paste0("TSNE_", i)]]
        names(obj) <- dat[["data.table"]][["cellNames"]]
        expect_true(identical(expected, obj))
    }
})

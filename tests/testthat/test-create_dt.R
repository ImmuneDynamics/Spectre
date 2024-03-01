#' Unit test for creating data.table from other popular SC based objects.
#' TODO FIX ME

# test_that("Non supported objects are not converted", {
#     to_convert <- c(1:20)
#     expect_error(create.dt(to_convert, from = "Vector"))
#     expect_error(create.dt(to_convert))
# })
# 
# test_that("Seurat RNA only object is converted", {
#     
#     pbmc <- qs::qread(test_path("testdata", "pbmc3k_seurat.qs"))
#     
# 
#     # Should succeed
#     expect_no_error(create.dt(pbmc, from = "Seurat"))
# 
#     dat <- create.dt(pbmc)
# 
#     # Meant to return a list
#     expect_equal(typeof(dat), "list")
# 
#     # Check all the genes are captured
#     # There is an element in the list that stores the gene names
#     # Check all the genes in each layer are captured
#     layers <- SeuratObject::Layers(pbmc)
#     
#     for (layer in layers) {
#         gene_names <- rownames(SeuratObject::LayerData(pbmc, assay = "RNA", layer=layer))
#         expect_true(all(gene_names %in% dat$featureNames$RNA[[layer]]))
#     }
#     
#     gene_names_in_cnt_mat <- lapply(layers, function(layer) {
#         gene_names <- rownames(SeuratObject::LayerData(pbmc, assay = "RNA", layer=layer))
#         paste(gene_names, "RNA", layer, sep = "_")
#     })
#     gene_names_in_cnt_mat <- unlist(gene_names_in_cnt_mat, use.names = FALSE)
#     expect_true(all(gene_names_in_cnt_mat %in% colnames(dat[["data.table"]])))
#     
#     
#     # Check all cell barcodes are captured
#     expect_true(all(colnames(pbmc) %in% dat[["data.table"]][["cellNames"]]))
#     expect_true(all(colnames(pbmc) %in% dat[["cellNames"]]))
# 
#     # Check metadata
#     expect_true(all(colnames(pbmc@meta.data) %in% dat[["meta.data"]]))
#     for (meta_dat in colnames(pbmc@meta.data)) {
#         expected <- pbmc@meta.data[[meta_dat]]
#         obj <- dat[["data.table"]][[meta_dat]]
#         expect_true(all.equal(obj, expected))
#     }
# 
#     # Check HVG
#     expect_true(all(Seurat::VariableFeatures(pbmc) %in% dat[["var.features"]]))
#     expect_true(all(head(Seurat::VariableFeatures(pbmc), 10) %in% dat[["var.features.top10"]]))
# 
#     # Check assay
#     expect_true(all(c("_RNA") %in% dat[["assays"]]))
# 
#     # Check all slots are captured
#     expect_true(all(c("counts", "data", "scale.data") %in% dat[["layers"]][['RNA']]))
# 
#     # Check dimreds
#     # TODO: find datasets that have dimreds
#     expect_true(all(c("pca_", "umap_") %in% dat[["dim.reds"]]))
#     # By default the PBMC data has 50 PCs.
#     pca_coordinates <- SeuratObject::Embeddings(pbmc, reduction = "pca")
#     for (i in c(1:50)) {
#         expected <- pca_coordinates[, paste0("PC_", i)]
#         obj <- dat[["data.table"]][[paste0("pca_PC_", i)]]
#         names(obj) <- dat[["data.table"]][["cellNames"]]
#         expect_true(identical(expected, obj))
#     }
# 
#     # By default, the PBMC data has 2 UMAP
#     umap_coordinates <- SeuratObject::Embeddings(pbmc, reduction = "umap")
#     for (i in c(1:2)) {
#         expected <- umap_coordinates[, paste0("umap_", i)]
#         obj <- dat[["data.table"]][[paste0("umap_umap_", i)]]
#         names(obj) <- dat[["data.table"]][["cellNames"]]
#         expect_true(identical(expected, obj))
#     }
# })
# 
# test_that("Seurat RNA+ADT only object is converted", {
#     
#     pbmc <- qs::qread(test_path("testdata", "cbmc_seurat.qs"))
#     
#     # Should succeed
#     expect_no_error(create.dt(pbmc, from = "Seurat"))
#     
#     dat <- create.dt(pbmc)
#     
#     # Meant to return a list
#     expect_equal(typeof(dat), "list")
#     
#     # Check all the genes are captured
#     # There is an element in the list that stores the gene names
#     # Check all the genes in each layer are captured
#     assays <- SeuratObject::Assays(pbmc)
#     
#     for (assay in assays) {
#         layers <- SeuratObject::Layers(pbmc[[assay]])
#         for (layer in layers) {
#             gene_names <- rownames(SeuratObject::LayerData(pbmc, assay = assay, layer=layer))
#             expect_true(all(gene_names %in% dat$featureNames[[assay]][[layer]]))
#         }
#     }
#     
#     
#     gene_names_in_cnt_mat <- lapply(assays, function(assay) {
#         layers <- SeuratObject::Layers(pbmc[[assay]])
#         gene_names_assay <- lapply(layers, function(layer) {
#             gene_names <- rownames(SeuratObject::LayerData(pbmc, assay = assay, layer=layer))
#             paste(gene_names, assay, layer, sep = "_")
#         })
#         unlist(gene_names_assay, use.names = FALSE)
#     })
#     gene_names_in_cnt_mat <- unlist(gene_names_in_cnt_mat, use.names = FALSE)
#     expect_true(all(gene_names_in_cnt_mat %in% colnames(dat[["data.table"]])))
#     
#     
#     # Check all cell barcodes are captured
#     expect_true(all(colnames(pbmc) %in% dat[["data.table"]][["cellNames"]]))
#     expect_true(all(colnames(pbmc) %in% dat[["cellNames"]]))
#     
#     # Check metadata
#     expect_true(all(colnames(pbmc@meta.data) %in% dat[["meta.data"]]))
#     for (meta_dat in colnames(pbmc@meta.data)) {
#         expected <- pbmc@meta.data[[meta_dat]]
#         obj <- dat[["data.table"]][[meta_dat]]
#         expect_true(all.equal(obj, expected))
#     }
#     
#     # Check HVG
#     expect_true(all(Seurat::VariableFeatures(pbmc) %in% dat[["var.features"]]))
#     expect_true(all(head(Seurat::VariableFeatures(pbmc), 10) %in% dat[["var.features.top10"]]))
#     
#     # Check assay
#     expect_true(all(c("_RNA") %in% dat[["assays"]]))
#     
#     # Check all slots are captured
#     expect_true(all(c("counts", "data", "scale.data") %in% dat[["layers"]][['RNA']]))
#     
#     # Check dimreds
#     expect_true(all(c("pca_", "umap_") %in% dat[["dim.reds"]]))
#     # By default the PBMC data has 50 PCs.
#     pca_coordinates <- SeuratObject::Embeddings(pbmc, reduction = "pca")
#     for (i in c(1:50)) {
#         expected <- pca_coordinates[, paste0("PC_", i)]
#         obj <- dat[["data.table"]][[paste0("pca_PC_", i)]]
#         names(obj) <- dat[["data.table"]][["cellNames"]]
#         expect_true(identical(expected, obj))
#     }
#     
#     # By default, the PBMC data has 2 UMAP
#     umap_coordinates <- SeuratObject::Embeddings(pbmc, reduction = "umap")
#     for (i in c(1:2)) {
#         expected <- umap_coordinates[, paste0("umap_", i)]
#         obj <- dat[["data.table"]][[paste0("umap_umap_", i)]]
#         names(obj) <- dat[["data.table"]][["cellNames"]]
#         expect_true(identical(expected, obj))
#     }
# })
# 
# 
# 
# test_that("SCE object is converted", {
#     
#     sce <- qs::qread(test_path("testdata", "sce_allendata.qs"))
# 
#     # Should succeed
#     expect_no_error(create.dt(sce, from = "SingleCellExperiment"))
# 
#     dat <- create.dt(sce)
# 
#     # Meant to return a list
#     expect_equal(typeof(dat), "list")
# 
#     # Check number of cells
#     expect_equal(dim(SummarizedExperiment::assay(sce))[2], nrow(dat[["data.table"]]))
# 
#     # Check all the genes are captured
#     gene_names <- c(
#         paste0(rownames(sce), "_counts"),
#         paste0(rownames(sce), "_logcounts")
#     )
#     expect_true(all(gene_names %in% colnames(dat[["data.table"]])))
#     # There is an element in the list that stores the gene names
#     expect_true(all(rownames(sce) %in% dat[["geneNames"]]))
# 
#     # Check all cell barcodes are captured
#     expect_true(all(colnames(sce) %in% dat[["data.table"]][["cellNames"]]))
#     expect_true(all(colnames(sce) %in% dat[["cellNames"]]))
# 
#     # Check metadata
#     expect_true(all(colnames(colData(sce)) %in% dat[["meta.data"]]))
#     for (meta_dat in colnames(colData(sce))) {
#         expected <- colData(sce)[[meta_dat]]
#         obj <- dat[["data.table"]][[meta_dat]]
#         expect_true(all.equal(obj, expected))
#     }
# 
#     # Check all counts assays are captured
#     expect_true(all(c("_counts", "_logcounts") %in% dat[["assays"]]))
# 
#     # Check dimreds
#     expect_true(all(c("PCA_", "TSNE_") %in% dat[["dim.reds"]]))
# 
#     # We set up the PCA to have 50 PCs
#     pca_coordinates <- reducedDim(sce, "PCA")
#     for (i in c(1:50)) {
#         expected <- pca_coordinates[, paste0("PC", i)]
#         obj <- dat[["data.table"]][[paste0("PCA_", i)]]
#         names(obj) <- dat[["data.table"]][["cellNames"]]
#         expect_true(identical(expected, obj))
#     }
# 
#     # By default, the data has 2 TSNE coordinates
#     tsne_coordinates <- reducedDim(sce, "TSNE")
#     for (i in c(1:2)) {
#         expected <- tsne_coordinates[, i]
#         obj <- dat[["data.table"]][[paste0("TSNE_", i)]]
#         names(obj) <- dat[["data.table"]][["cellNames"]]
#         expect_true(identical(expected, obj))
#     }
# })

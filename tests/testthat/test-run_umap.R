test_that("run.umap works with fast = TRUE", {

    dat <- do.subsample(demo.clustered, 100)

    use.cols <- c(
        "NK11_asinh", "CD3_asinh", "CD45_asinh", "Ly6G_asinh", 
        "CD11b_asinh", "B220_asinh", "CD8a_asinh", "Ly6C_asinh", "CD4_asinh"
    )

    result <- run.umap(
        dat = dat, 
        use.cols = use.cols, 
        fast = TRUE,
        umap.x.name = "UMAP_1",
        umap.y.name = "UMAP_2",
        verbose = FALSE
    )

    expect_s3_class(result, "data.table")
    expect_true("UMAP_1" %in% names(result))
    expect_true("UMAP_2" %in% names(result))
    expect_equal(nrow(result), nrow(dat))

    expect_false("UMAP_1" %in% names(dat))
    expect_false("UMAP_2" %in% names(dat))
})


test_that("run.umap works with fast = FALSE (umap)", {
    dat <- do.subsample(demo.clustered, 100)

    use.cols <- c(
        "NK11_asinh", "CD3_asinh", "CD45_asinh", "Ly6G_asinh", 
        "CD11b_asinh", "B220_asinh", "CD8a_asinh", "Ly6C_asinh", "CD4_asinh"
    )

    result <- run.umap(
        dat = dat, 
        use.cols = use.cols, 
        fast = FALSE,
        umap.x.name = "UMAP_1",
        umap.y.name = "UMAP_2",
        verbose = FALSE
    )

    expect_s3_class(result, "data.table")
    expect_true("UMAP_1" %in% names(result))
    expect_true("UMAP_2" %in% names(result))
    expect_equal(nrow(result), nrow(dat))

    expect_false("UMAP_1" %in% names(dat))
    expect_false("UMAP_2" %in% names(dat))
})

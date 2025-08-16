test_that("run.flowsom() runs on valid input", {
    dt <- demo.clustered

    use.cols <- c("NK11_asinh", "CD3_asinh", "CD45_asinh")

    result <- run.flowsom(
        dat = dt,
        use.cols = use.cols,
        xdim = 5,
        ydim = 5,
        clust.seed = 123,
        meta.seed = 321,
        clust.name = "Cluster",
        meta.clust.name = "Meta",
        meta.k = 5,
        verbose = FALSE
    )

    # Check that the cluster columns exist
    expect_true("Cluster" %in% names(result))
    expect_true("Meta" %in% names(result))

    # dt shouldn't have changed
    expect_false("Cluster" %in% names(dt))
    expect_false("Meta" %in% names(dt))

    # Check that both are integer vectors of correct length
    expect_equal(length(result$Cluster), nrow(dt))
    expect_equal(length(result$Meta), nrow(dt))

    # Check original columns are still present
    expect_true(all(names(dt) %in% names(result)))
})

test_that("run.flowsom() returns reproducible clustering", {
    dt <- demo.clustered

    use.cols <- c("NK11_asinh", "CD3_asinh", "CD45_asinh")

    res1 <- run.flowsom(
        dat = dt,
        use.cols = use.cols,
        clust.seed = 42,
        meta.seed = 99,
        meta.k = 5,
        clust.name = "Cluster",
        meta.clust.name = "Meta",
        verbose = FALSE
    )

    res2 <- run.flowsom(
        dat = dt,
        use.cols = use.cols,
        clust.seed = 42,
        meta.seed = 99,
        meta.k = 5,
        clust.name = "Cluster",
        meta.clust.name = "Meta",
        verbose = FALSE
    )

    expect_equal(res1$Cluster, res2$Cluster)
    expect_equal(res1$Meta, res2$Meta)

    # dt shouldn't have changed
    expect_false("Cluster" %in% names(dt))
    expect_false("Meta" %in% names(dt))
})

test_that("run.flowsom() skips metaclustering when meta.k = 0", {
    dt <- demo.clustered

    use.cols <- c("NK11_asinh", "CD3_asinh", "CD45_asinh")

    result <- run.flowsom(
        dat = dt,
        use.cols = use.cols,
        meta.k = 0,
        clust.name = "Cluster",
        meta.clust.name = "Meta",
        verbose = FALSE
    )

    expect_true("Cluster" %in% names(result))
    expect_false("Meta" %in% names(result))
})

test_that("run.flowsom() correctness", {
    dt <- demo.clustered

    use.cols <- c(
        "NK11_asinh", "CD3_asinh", "CD45_asinh", 
        "Ly6G_asinh", "CD11b_asinh", "B220_asinh", 
        "CD8a_asinh", "Ly6C_asinh", "CD4_asinh"
    )

    # same seeting as how demo clustered data was generated
    result <- run.flowsom(
        dat = dt,
        use.cols = use.cols,
        xdim = 14,
        ydim = 14,
        meta.k = 7,
        clust.seed = 42,
        meta.seed = 42,
        clust.name = "Cluster",
        meta.clust.name = "Meta",
        verbose = FALSE
    )

    expect_equal(as.numeric(result$Cluster), as.numeric(dt$FlowSOM_cluster))
    expect_equal(as.numeric(result$Meta), as.numeric(dt$FlowSOM_metacluster))

})



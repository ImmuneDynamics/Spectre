# Use the demo data. Subsampled so it is quicker to run.

test_that("Fast UMAP single thread can run", {
    dat <- do.subsample(
        dat = demo.data, 
        targets = 1000, 
        seed = 42
    )
    markers <- names(dat)[11:19]
    
    dat <- run_actual_umap(
        dat, 
        use.cols = markers, 
        umap.x.name = "UWOT_X", 
        umap.y.name = "UWOT_Y",
        fast = TRUE,
        verbose = FALSE,
        n_threads = 1
    )
    
    expect_true("UWOT_X" %in% names(dat))
    expect_true("UWOT_Y" %in% names(dat))

    expect_equal(
        nrow(na.omit(dat, cols = c("UWOT_X", "UWOT_Y"))),
        nrow(dat)
    )
})

test_that("Fast UMAP infer n_threads can run", {
    dat <- do.subsample(
        dat = demo.data, 
        targets = 1000, 
        seed = 42
    )
    markers <- names(dat)[11:19]
    
    dat <- run_actual_umap(
        dat, 
        use.cols = markers, 
        umap.x.name = "UWOT_X", 
        umap.y.name = "UWOT_Y",
        fast = TRUE,
        verbose = FALSE,
        n_threads = 'auto'
    )
    
    expect_true("UWOT_X" %in% names(dat))
    expect_true("UWOT_Y" %in% names(dat))
    
    expect_equal(
        nrow(na.omit(dat, cols = c("UWOT_X", "UWOT_Y"))),
        nrow(dat)
    )
})

test_that("Default UMAP can run", {
    dat <- do.subsample(
        dat = demo.data, 
        targets = 1000, 
        seed = 42
    )
    markers <- names(dat)[11:19]
    
    dat <- run_actual_umap(
        dat, 
        use.cols = markers, 
        umap.x.name = "UMAP_X_2", 
        umap.y.name = "UMAP_Y_2",
        fast = FALSE,
        verbose = FALSE
    )
    
    expect_true("UMAP_X_2" %in% names(dat))
    expect_true("UMAP_Y_2" %in% names(dat))
    
    expect_equal(
        nrow(na.omit(dat, cols = c("UMAP_X_2", "UMAP_Y_2"))),
        nrow(dat)
    )
})

# ---- run.umap test ----

test_that("run.umap on data.table can run", {
    dat <- do.subsample(
        dat = demo.data, 
        targets = 1000, 
        seed = 42
    )
    markers <- names(dat)[11:19]
    
    dat <- run.umap(
        dat, 
        use.cols = markers, 
        umap.x.name = "UMAP_X_2", 
        umap.y.name = "UMAP_Y_2",
        fast = FALSE,
        verbose = FALSE
    )
    
    expect_true("UMAP_X_2" %in% names(dat))
    expect_true("UMAP_Y_2" %in% names(dat))
    
})

test_that("run.umap fast on data.table can run", {
    dat <- do.subsample(
        dat = demo.data, 
        targets = 1000, 
        seed = 42
    )
    markers <- names(dat)[11:19]
    
    dat <- run.umap(
        dat, 
        use.cols = markers, 
        umap.x.name = "UMAP_X_FAST", 
        umap.y.name = "UMAP_Y_FAST",
        fast = TRUE,
        verbose = FALSE
    )
    
    expect_true("UMAP_X_FAST" %in% names(dat))
    expect_true("UMAP_Y_FAST" %in% names(dat))
})


test_that("run.umap on Spectre object can run", {
    dat <- do.subsample(
        dat = demo.data, 
        targets = 1000, 
        seed = 42
    )
    markers <- names(dat)[11:19]
    dat[, cell_id := paste0("Cell_", seq(1000))]
    
    obj <- create.spectre.object(cell_id_col = "cell_id")
    obj <- add.new.data(obj, dat, "test")
    
    obj <- run.umap(
        obj,
        data_source = "test",
        output_name = "test_umap",
        use.cols = markers,
        umap.x.name = "UMAP_X_2",
        umap.y.name = "UMAP_Y_2",
        fast = FALSE,
        verbose = FALSE
    )
    
    expect_true("test_umap" %in% names(obj))
    expect_true("UMAP_X_2" %in% names(obj$test_umap))
    expect_true("UMAP_Y_2" %in% names(obj$test_umap))
    expect_equal(ncol(obj$test_umap), 3)
    # content should not be modified
    expect_equal(ncol(obj$test), ncol(dat))
})

test_that("run.umap fast on Spectre object can run", {
    dat <- do.subsample(
        dat = demo.data, 
        targets = 1000, 
        seed = 42
    )
    markers <- names(dat)[11:19]
    dat[, cell_id := paste0("Cell_", seq(1000))]
    
    obj <- create.spectre.object(cell_id_col = "cell_id")
    obj <- add.new.data(obj, dat, "test")
    
    obj <- run.umap(
        obj,
        data_source = "test",
        output_name = "test_umap",
        use.cols = markers,
        umap.x.name = "UMAP_X_2",
        umap.y.name = "UMAP_Y_2",
        fast = TRUE,
        verbose = FALSE
    )
    
    expect_true("test_umap" %in% names(obj))
    expect_true("UMAP_X_2" %in% names(obj$test_umap))
    expect_true("UMAP_Y_2" %in% names(obj$test_umap))
    expect_equal(ncol(obj$test_umap), 3)
    # content should not be modified
    expect_equal(ncol(obj$test), ncol(dat))
})

test_that("run.umap update data_source on Spectre object can run", {
    dat <- do.subsample(
        dat = demo.data, 
        targets = 1000, 
        seed = 42
    )
    markers <- names(dat)[11:19]
    dat[, cell_id := paste0("Cell_", seq(1000))]
    
    obj <- create.spectre.object(cell_id_col = "cell_id")
    obj <- add.new.data(obj, dat, "test")
    
    obj_umap <- run.umap(
        obj,
        data_source = "test",
        output_name = NULL,
        use.cols = markers,
        umap.x.name = "UMAP_X_2",
        umap.y.name = "UMAP_Y_2",
        fast = FALSE,
        verbose = FALSE
    )
    
    expect_true("UMAP_X_2" %in% names(obj_umap$test))
    expect_true("UMAP_Y_2" %in% names(obj_umap$test))
    expect_equal(ncol(obj_umap$test), ncol(obj$test) + 2)
})

test_that("run.umap fast update data_source on Spectre object can run", {
    dat <- do.subsample(
        dat = demo.data, 
        targets = 1000, 
        seed = 42
    )
    markers <- names(dat)[11:19]
    dat[, cell_id := paste0("Cell_", seq(1000))]
    
    obj <- create.spectre.object(cell_id_col = "cell_id")
    obj <- add.new.data(obj, dat, "test")
    
    obj_umap <- run.umap(
        obj,
        data_source = "test",
        output_name = NULL,
        use.cols = markers,
        umap.x.name = "UMAP_X_2",
        umap.y.name = "UMAP_Y_2",
        fast = TRUE,
        verbose = FALSE
    )
    
    expect_true("UMAP_X_2" %in% names(obj_umap$test))
    expect_true("UMAP_Y_2" %in% names(obj_umap$test))
    expect_equal(ncol(obj_umap$test), ncol(obj$test) + 2)
})

test_that("run.umap same column names are overwritten", {
    dat <- do.subsample(
        dat = demo.data, 
        targets = 1000, 
        seed = 42
    )
    markers <- names(dat)[11:19]
    dat[, cell_id := paste0("Cell_", seq(1000))]
    
    dat_umap <- run.umap(
        dat, 
        use.cols = markers, 
        umap.x.name = "UMAP_X", 
        umap.y.name = "UMAP_Y",
        fast = TRUE,
        verbose = FALSE
    )
    
    expect_false(setequal(dat$UMAP_X, dat_umap$UMAP_X))
    expect_false(setequal(dat$UMAP_Y, dat_umap$UMAP_Y))
    
    obj <- create.spectre.object(cell_id_col = "cell_id")
    obj <- add.new.data(obj, dat, "test")
    
    obj_umap <- run.umap(
        obj,
        data_source = "test",
        output_name = NULL,
        use.cols = markers,
        umap.x.name = "UMAP_X",
        umap.y.name = "UMAP_Y",
        fast = TRUE,
        verbose = FALSE
    )
    
    expect_false(setequal(obj$test$UMAP_X, obj_umap$test$UMAP_X))
    expect_false(setequal(obj$test$UMAP_Y, obj_umap$test$UMAP_Y))
})



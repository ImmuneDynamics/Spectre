# All do_actual_transformation tests that should work

test_that("do_actual_transformation single cofactor works correctly", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    cofactor <- 7
    
    res <- do_actual_transformation(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = cofactor,
        append.cf = FALSE,
        verbose = FALSE
    )
    actual_dat <- res$transformed_val
    
    expected_dat <- asinh(dat / cofactor)
    
    expect_equal(actual_dat$marker1_asinh, expected_dat$marker1)
    expect_equal(actual_dat$marker2_asinh, expected_dat$marker2)
    
    
})


test_that("do_actual_transformation multiple cofactors works correctly", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    cofactor <- c(5,10)
    
    res <- do_actual_transformation(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = cofactor,
        append.cf = FALSE,
        verbose = FALSE
    )
    actual_dat <- res$transformed_val
    
    expected_dat <- data.table(
        marker1_asinh = asinh(dat$marker1 / 5),
        marker2_asinh = asinh(dat$marker2 / 10)
    )
    
    expect_equal(actual_dat$marker1_asinh, expected_dat$marker1_asinh)
    expect_equal(actual_dat$marker2_asinh, expected_dat$marker2_asinh)
    
    
})


test_that("do_actual_transformation use.cols works correctly", {
    # marker3 should not be transformed
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2),
        marker3 = rnorm(10, 5)
    )
    
    res <- do_actual_transformation(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = 10,
        append.cf = FALSE,
        verbose = FALSE
    )
    
    expect_equal(names(res$transformed_val), c("marker1_asinh", "marker2_asinh"))
    
})

test_that("do_actual_transformation append_cf = TRUE works correctly", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    
    res <- do_actual_transformation(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = 7, 
        append.cf = TRUE,
        verbose = FALSE
    )
    
    expect_equal(names(res$transformed_val), c("marker1_asinh_cf7", "marker2_asinh_cf7"))
    
    res <- do_actual_transformation(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = c(5,10), 
        append.cf = TRUE,
        verbose = FALSE
    )
    
    expect_equal(names(res$transformed_val), c("marker1_asinh_cf5", "marker2_asinh_cf10"))
    
})

test_that("do_actual_transformation append_cf = FALSE works correctly", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    
    # single cofactor
    res <- do_actual_transformation(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = 7, 
        append.cf = FALSE,
        verbose = FALSE
    )
    
    expect_equal(names(res$transformed_val), c("marker1_asinh", "marker2_asinh"))
    
    # multiple cofactor
    res <- do_actual_transformation(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = c(5,10), 
        append.cf = FALSE,
        verbose = FALSE
    )
    
    expect_equal(names(res$transformed_val), c("marker1_asinh", "marker2_asinh"))
    
    # automated cofactor
    res <- do_actual_transformation(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = NULL,
        cofactor_inference_method = "flowVS", 
        append.cf = FALSE,
        verbose = FALSE
    )
    
    expect_equal(names(res$transformed_val), c("marker1_asinh", "marker2_asinh"))
    
    res <- do_actual_transformation(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = NULL,
        cofactor_inference_method = "top10",
        append.cf = FALSE,
        verbose = FALSE
    )
    
    expect_equal(names(res$transformed_val), c("marker1_asinh", "marker2_asinh"))
})

test_that("do_actual_transformation round = TRUE works correctly", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    
    cofactor <- 7
    
    res <- do_actual_transformation(
        dat = dat, 
        use.cols = c("marker1", "marker2"), 
        cofactor = cofactor,
        append.cf = FALSE, 
        digits = 3,
        verbose = FALSE
    )
    actual_dat <- res$transformed_val
    
    expect_equal(names(actual_dat), c("marker1_asinh", "marker2_asinh"))
    
    expected_dat <- round(asinh(dat / cofactor), 3)
    
    expect_equal(actual_dat$marker1_asinh, expected_dat$marker1)
    expect_equal(actual_dat$marker2_asinh, expected_dat$marker2)
    
    # for automated cofactor calculation
    expect_no_error(
        do_actual_transformation(
            dat = dat, 
            use.cols = c("marker1", "marker2"), 
            cofactor = NULL,
            cofactor_inference_method = "flowVS",
            append.cf = FALSE, 
            digits = 3,
            verbose = FALSE
        )
    )
    
    expect_no_error(
        do_actual_transformation(
            dat = dat, 
            use.cols = c("marker1", "marker2"), 
            cofactor = NULL,
            cofactor_inference_method = "top10",
            append.cf = FALSE, 
            digits = 3,
            verbose = FALSE
        )
    )
    
    
})


# All do.asinh tests that should succeed

test_that("do.asinh using flowVS works", {
    dat <- data.table(
        cell_id = paste0("cell_", seq(1, 10)),
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2),
        marker3 = rnorm(10, 3)
    )
    
    dat_actual = do.asinh(
        dat = dat, 
        use.cols = c("marker1", "marker2", "marker3"), 
        cofactor = NULL,
        cofactor_inference_method = "flowVS",
        append.cf = FALSE,
        verbose = FALSE
    )
    
    expect_equal(names(dat_actual), c("cell_id", "marker1", "marker2", "marker3", 
                                      "marker1_asinh", "marker2_asinh", "marker3_asinh"))
    expect_equal(nrow(dat_actual), 10)
    
    # for spectre object
    obj <- create.spectre.object(cell_id_col = "cell_id")
    obj <- add.new.data(obj, dat, "test")
    
    obj_actual = do.asinh(
        dat = obj,
        data_source = "test",
        output_name = "test_auto_asinh",
        use.cols = c("marker1", "marker2", "marker3"), 
        cofactor = NULL,
        cofactor_inference_method = "flowVS",
        verbose = FALSE,
        append.cf = FALSE
    )
    
    # check the actual raw data is not touched!
    expect_equal(names(obj_actual$test), c("cell_id","marker1", "marker2", "marker3"))
    expect_equal(names(obj_actual$test_auto_asinh), c("cell_id", "marker1", "marker2", "marker3"))
    expect_equal(nrow(obj_actual$test_auto_asinh), 10)
    expect_equal(nrow(
        attributes(obj_actual$test_auto_asinh)$cofactors
    ), 3)
    
    
})

test_that("do.asinh using top10 works", {
    dat <- data.table(
        cell_id = paste0("cell_", seq(1, 10)),
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2),
        marker3 = rnorm(10, 3)
    )
    
    dat_actual = do.asinh(
        dat = dat, 
        use.cols = c("marker1", "marker2", "marker3"), 
        cofactor = NULL,
        cofactor_inference_method = "top10",
        append.cf = FALSE,
        verbose = FALSE
    )
    
    expect_equal(names(dat_actual), c("cell_id", "marker1", "marker2", "marker3", 
                                      "marker1_asinh", "marker2_asinh", "marker3_asinh"))
    expect_equal(nrow(dat_actual), 10)
    
    # for spectre object
    obj <- create.spectre.object(cell_id_col = "cell_id")
    obj <- add.new.data(obj, dat, "test")
    
    obj_actual = do.asinh(
        dat = obj,
        data_source = "test",
        output_name = "test_auto_asinh",
        use.cols = c("marker1", "marker2", "marker3"), 
        cofactor = NULL,
        cofactor_inference_method = "top10",
        verbose = FALSE,
        append.cf = FALSE
    )
    
    # check the actual raw data is not touched!
    expect_equal(names(obj_actual$test), c("cell_id","marker1", "marker2", "marker3"))
    expect_equal(names(obj_actual$test_auto_asinh), c("cell_id", "marker1", "marker2", "marker3"))
    expect_equal(nrow(obj_actual$test_auto_asinh), 10)
    expect_equal(nrow(
        attributes(obj_actual$test_auto_asinh)$cofactors
    ), 3)
    
    
    
})

test_that("do.asinh single cofactor works", {
    dat <- data.table(
        cell_id = paste0("cell_", seq(1, 10)),
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    cofactor <- 7
    
    dat_actual <- do.asinh(dat, use.cols = c("marker1", "marker2"), 
                    cofactor = cofactor, verbose = FALSE,
                    append.cf = FALSE)
    
    expect_equal(names(dat_actual), c("cell_id", "marker1", "marker2", "marker1_asinh", "marker2_asinh"))
    expect_equal(nrow(dat_actual), 10)
    
    obj <- create.spectre.object(cell_id_col = "cell_id")
    obj <- add.new.data(obj, dat, "test")
    
    actual_obj <- do.asinh(
        dat = obj, 
        data_source = "test",
        output_name = "test_asinh",
        use.cols = c("marker1", "marker2"), 
        cofactor = 7, 
        verbose = FALSE,
        append.cf = FALSE
    )
    
    # check the actual raw data is not touched!
    expect_equal(names(actual_obj$test), c("cell_id","marker1", "marker2"))
    expect_identical(actual_obj$test, dat)
    expect_equal(names(actual_obj$test_asinh), c("cell_id", "marker1", "marker2"))
    expect_equal(nrow(actual_obj$test_asinh), 10)
})



# All do.asinh tests that meant to fail

test_that("do.asinh missing use.cols fails", {
    dat <- data.table(
        cell_id = paste0("cell_", seq(1, 10)),
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )

    expect_error(do.asinh(dat, verbose = FALSE))

    obj <- create.spectre.object(cell_id_col = "cell_id")
    obj <- add.new.data(obj, dat, "test")
    expect_error(do.asinh(obj, verbose = FALSE))
})


test_that("do.asinh non-numeric use.cols fails", {
    dat <- data.table(
        cell_id = paste0("cell_", seq(1, 10)),
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2),
        marker3 = rep("A", 10)
    )

    expect_error(do.asinh(dat, use.cols = c("marker1", "marker3") ,verbose = FALSE))

    obj <- create.spectre.object(cell_id_col = "cell_id")
    obj <- add.new.data(obj, dat, "test")
    expect_error(do.asinh(obj, use.cols = c("marker1", "marker3") ,verbose = FALSE))
})


test_that("do.asinh missing some use.cols columns fails", {
    dat <- data.table(
        cell_id = paste0("cell_", seq(1, 10)),
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )

    expect_error(do.asinh(dat, use.cols = c("marker1", "marker4") ,verbose = FALSE))

    obj <- create.spectre.object(cell_id_col = "cell_id")
    obj <- add.new.data(obj, dat, "test")
    expect_error(do.asinh(obj, use.cols = c("marker1", "marker4") ,verbose = FALSE))
})

test_that("do.asinh more use.cols than multiple cofactors errors", {
    dat <- data.table(
        cell_id = paste0("cell_", seq(1, 10)),
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2),
        marker3 = rnorm(10, 3)
    )
    
    expect_error(
        do.asinh(
            dat = dat, 
            use.cols = c("marker3"), 
            cofactor = c(5,10), 
            append.cf = TRUE,
            verbose = FALSE
        )
    )
    obj <- create.spectre.object(cell_id_col = "cell_id")
    obj <- add.new.data(obj, dat, "test")
    
    expect_error(
        do.asinh(
            dat = obj, 
            use.cols = c("marker3"), 
            cofactor = c(5,10), 
            append.cf = TRUE,
            verbose = FALSE
        )
    )
})

test_that("do.asinh non specified automated inferrence failed", {
    dat <- data.table(
        cell_id = paste0("cell_", seq(1, 10)),
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2),
        marker3 = rnorm(10, 3)
    )
    
    expect_error(
        do.asinh(
            dat = dat, 
            use.cols = c("marker3"), 
            cofactor = NULL,
            cofactor_inference_method = "hello",
            append.cf = TRUE,
            verbose = FALSE
        )
    )
    obj <- create.spectre.object(cell_id_col = "cell_id")
    obj <- add.new.data(obj, dat, "test")
    
    expect_error(
        do.asinh(
            dat = obj, 
            use.cols = c("marker3"), 
            cofactor = NULL,
            cofactor_inference_method = "hello",
            append.cf = TRUE,
            verbose = FALSE
        )
    )
})







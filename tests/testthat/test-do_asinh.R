test_that("do_actual_transformation single cofactor works correctly", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    cofactor <- 7
    
    actual_dat <- do_actual_transformation(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = cofactor,
        append.cf = FALSE,
        verbose = FALSE
    )
    
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
    
    actual_dat <- do_actual_transformation(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = cofactor,
        append.cf = FALSE,
        verbose = FALSE
    )
    
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
    
    actual_dat <- do_actual_transformation(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = 10,
        append.cf = FALSE,
        verbose = FALSE
    )
    
    expect_equal(names(actual_dat), c("marker1_asinh", "marker2_asinh"))
    
})

test_that("do_actual_transformation with append_cf works correctly", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    
    actual_dat <- do_actual_transformation(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = 7, 
        append.cf = TRUE,
        verbose = FALSE
    )
    
    expect_equal(names(actual_dat), c("marker1_asinh_cf7", "marker2_asinh_cf7"))
    
    actual_dat <- do_actual_transformation(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = c(5,10), 
        append.cf = TRUE,
        verbose = FALSE
    )
    
    expect_equal(names(actual_dat), c("marker1_asinh_cf5", "marker2_asinh_cf10"))
    
})

test_that("do_actual_transformation with round works correctly", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    
    cofactor <- 7
    
    actual_dat <- do_actual_transformation(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = cofactor,
        append.cf = FALSE, digits = 3,
        verbose = FALSE
    )
    
    expect_equal(names(actual_dat), c("marker1_asinh", "marker2_asinh"))
    
    expected_dat <- round(asinh(dat / cofactor), 3)
    
    expect_equal(actual_dat$marker1_asinh, expected_dat$marker1)
    expect_equal(actual_dat$marker2_asinh, expected_dat$marker2)
    
    
})

test_that("do_actual_transformation more use.cols than cofactor errors", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2),
        marker3 = rnorm(10, 3)
    )
    
    expect_error(
        do_actual_transformation(
            dat = dat, 
            use.cols = c("marker3"), 
            cofactor = c(5,10), 
            append.cf = TRUE,
            verbose = FALSE
        )
    )
    
    
})

test_that("do.asinh for data.table works", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    cofactor <- 7
    
    dat <- do.asinh(dat, use.cols = c("marker1", "marker2"), 
                    cofactor = cofactor, verbose = FALSE,
                    append.cf = FALSE)
    
    expect_equal(names(dat), c("marker1", "marker2", "marker1_asinh", "marker2_asinh"))
    expect_equal(nrow(dat), 10)
})

test_that("do.asinh for Spectre object works", {
    dat <- data.table(
        cell_id = paste0("cell_", seq(1, 10)),
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
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
    expect_equal(names(actual_obj$test_asinh), c("cell_id", "marker1_asinh", "marker2_asinh"))
    expect_equal(nrow(actual_obj$test_asinh), 10)
})

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


test_that("asinh non-numeric use.cols fails", {
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


test_that("asinh missing some use.cols columns fails", {
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











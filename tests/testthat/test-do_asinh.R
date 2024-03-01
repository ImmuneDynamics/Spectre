test_that("do_actual_transformation works correctly", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    cofactor <- 7
    
    actual_dat <- do_actual_transformation(dat, use.cols = c("marker1", "marker2"), 
                    cofactor = cofactor,
                    append.cf = FALSE)
    
    expected_dat <- asinh(dat / cofactor)
    
    expect_equal(actual_dat$marker1_asinh, expected_dat$marker1)
    expect_equal(actual_dat$marker2_asinh, expected_dat$marker2)
    
    
})


test_that("do_actual_transformation use.cols works correctly", {
    # marker3 should not be transformed
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2),
        marker3 = rnorm(10, 5)
    )
    
    actual_dat <- do_actual_transformation(dat, use.cols = c("marker1", "marker2"), 
                    cofactor = 10,
                    append.cf = FALSE)
    
    expect_equal(names(actual_dat), c("marker1_asinh", "marker2_asinh"))
    
})

test_that("do_actual_transformation with append_cf works correctly", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    
    actual_dat <- do_actual_transformation(dat, use.cols = c("marker1", "marker2"), 
                    cofactor = 7, 
                    append.cf = TRUE)
    
    expect_equal(names(actual_dat), c("marker1_asinh_cf7", "marker2_asinh_cf7"))
    
})

test_that("do_actual_transformation with round works correctly", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    
    cofactor <- 7
    
    actual_dat <- do_actual_transformation(dat, use.cols = c("marker1", "marker2"), 
                    cofactor = cofactor,
                    append.cf = FALSE, digits = 3)
    
    expect_equal(names(actual_dat), c("marker1_asinh", "marker2_asinh"))
    
    expected_dat <- round(asinh(dat / cofactor), 3)
    
    expect_equal(actual_dat$marker1_asinh, expected_dat$marker1)
    expect_equal(actual_dat$marker2_asinh, expected_dat$marker2)
    
    
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

test_that("do.asinh for SpectreObject works", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    obj <- SpectreObject(cytometry_data=dat)
    
    actual_obj <- do.asinh(obj, slot_name = "cytometry_data",
                    use.cols = c("marker1", "marker2"), 
                    cofactor = 7, verbose = FALSE,
                    append.cf = FALSE)
    
    expect_equal(names(actual_obj@cytometry_data), c("marker1", "marker2", "marker1_asinh", "marker2_asinh"))
    expect_equal(nrow(actual_obj@cytometry_data), 10)
})

test_that("do.asinh missing use.cols fails", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    
    expect_error(do.asinh(dat, verbose = FALSE))
    
    obj <- SpectreObject(cytometry_data=dat)
    expect_error(do.asinh(obj, verbose = FALSE))
})


test_that("asinh no numeric use.cols fails", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2),
        marker3 = rep("A", 10)
    )
    
    expect_error(do.asinh(dat, use.cols = c("marker1", "marker3") ,verbose = FALSE))
    
    obj <- SpectreObject(cytometry_data=dat)
    expect_error(do.asinh(obj, use.cols = c("marker1", "marker3") ,verbose = FALSE))
})


test_that("asinh missing some use.cols columns fails", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    
    expect_error(do.asinh(dat, use.cols = c("marker1", "marker4") ,verbose = FALSE))
    
    obj <- SpectreObject(cytometry_data=dat)
    expect_error(do.asinh(obj, use.cols = c("marker1", "marker4") ,verbose = FALSE))
})











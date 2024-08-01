# All do.asinh tests that should work

test_that("do.asinh single cofactor works correctly", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    cofactor <- 7
    
    res <- do.asinh(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = cofactor,
        append.cf = FALSE,
        verbose = FALSE
    )
    expected_dat <- asinh(dat / cofactor)
    
    expect_equal(res$marker1_asinh, expected_dat$marker1)
    expect_equal(res$marker2_asinh, expected_dat$marker2)
    
    # check cofactors are appended
    expect_equal(attributes(res)$cofactors[['marker1']], 7)
    expect_equal(attributes(res)$cofactors[['marker2']], 7)
    
    
})


test_that("do.asinh multiple cofactors works correctly", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    cofactor <- c(5,10)
    
    res <- do.asinh(
        dat, use.cols = c("marker1", "marker2"), 
        cofactor = cofactor,
        append.cf = FALSE,
        verbose = FALSE
    )

    expected_dat <- data.table(
        marker1_asinh = asinh(dat$marker1 / 5),
        marker2_asinh = asinh(dat$marker2 / 10)
    )
    
    expect_equal(res$marker1_asinh, expected_dat$marker1_asinh)
    expect_equal(res$marker2_asinh, expected_dat$marker2_asinh)
    
    # check cofactors are appended
    expect_equal(attributes(res)$cofactors[['marker1']], 5)
    expect_equal(attributes(res)$cofactors[['marker2']], 10)
    
    
})


test_that("do.asinh use.cols works correctly", {
    # marker3 should not be transformed
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2),
        marker3 = rnorm(10, 5)
    )

    res <- do.asinh(
        dat, use.cols = c("marker1", "marker2"),
        cofactor = 10,
        append.cf = FALSE,
        verbose = FALSE,
        add_to_table = FALSE
    )

    expect_equal(names(res), c("marker1_asinh", "marker2_asinh"))

})

test_that("do.asinh append_cf = TRUE works correctly", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )

    res <- do.asinh(
        dat, use.cols = c("marker1", "marker2"),
        cofactor = 7,
        append.cf = TRUE,
        verbose = FALSE,
        add_to_table = FALSE
    )

    expect_equal(names(res), c("marker1_asinh_cf7", "marker2_asinh_cf7"))

    res <- do.asinh(
        dat, use.cols = c("marker1", "marker2"),
        cofactor = c(5,10),
        append.cf = TRUE,
        verbose = FALSE,
        add_to_table = TRUE
    )

    expect_true(all(c("marker1_asinh_cf5", "marker2_asinh_cf10") %in% names(res)))

})

test_that("do.asinh round = TRUE works correctly", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )

    cofactor <- 7

    res <- do.asinh(
        dat = dat,
        use.cols = c("marker1", "marker2"),
        cofactor = cofactor,
        append.cf = FALSE,
        digits = 3,
        verbose = FALSE,
        add_to_table = FALSE
    )

    expect_equal(names(res), c("marker1_asinh", "marker2_asinh"))

    expected_dat <- round(asinh(dat / cofactor), 3)

    expect_equal(res$marker1_asinh, expected_dat$marker1)
    expect_equal(res$marker2_asinh, expected_dat$marker2)

})


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

    expect_equal(length(names(dat_actual)), 7)
    expect_equal(nrow(dat_actual), 10)
    
    # check cofactors are appended as attributes
    expect_no_error(attributes(dat_actual)$marker1)
    expect_no_error(attributes(dat_actual)$marker2)
    expect_no_error(attributes(dat_actual)$marker3)

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

    expect_equal(length(names(dat_actual)), 7)
    expect_equal(nrow(dat_actual), 10)
    
    # check cofactors are appended as attributes
    expect_no_error(attributes(dat_actual)$marker1)
    expect_no_error(attributes(dat_actual)$marker2)
    expect_no_error(attributes(dat_actual)$marker3)

})

# All do.asinh tests that meant to fail

test_that("do.asinh missing use.cols fails", {
    dat <- data.table(
        cell_id = paste0("cell_", seq(1, 10)),
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )

    expect_error(do.asinh(dat, verbose = FALSE))

})


test_that("do.asinh non-numeric use.cols fails", {
    dat <- data.table(
        cell_id = paste0("cell_", seq(1, 10)),
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2),
        marker3 = rep("A", 10)
    )

    expect_error(do.asinh(dat, use.cols = c("marker1", "marker3") ,verbose = FALSE))

})


test_that("do.asinh missing some use.cols columns fails", {
    dat <- data.table(
        cell_id = paste0("cell_", seq(1, 10)),
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )

    expect_error(do.asinh(dat, use.cols = c("marker1", "marker4") ,verbose = FALSE))

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
})


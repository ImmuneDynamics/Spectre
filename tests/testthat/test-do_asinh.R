test_that("do.asinh works with basic input", {
    # Sample data
    dat <- data.table(
        A = c(0, 5, 10),
        B = c(1, 2, 3),
        C = c("x", "y", "z")  # Non-numeric column
    )
    
    # Apply transformation
    result <- do.asinh(dat, use.cols = c("A", "B"))

    # Check that original columns are still there
    expect_true(all(c("A", "B", "C") %in% names(result)))

    # Check that transformed columns are added
    expect_true(all(c("A_asinh", "B_asinh") %in% names(result)))

    # Manually calculate expected transformation
    expected_A <- asinh(dat$A / 5)
    expected_B <- asinh(dat$B / 5)

    # Compare results
    expect_equal(result$A_asinh, expected_A)
    expect_equal(result$B_asinh, expected_B)
})

test_that("do.asinh appends cofactor to new column names when append.cf = TRUE", {
    dat <- data.table(A = c(1, 2, 3))
    result <- do.asinh(dat, use.cols = "A", cofactor = 10, append.cf = TRUE)

    expect_true("A_asinh_cf10" %in% names(result))
    expect_equal(result$A_asinh_cf10, asinh(dat$A / 10))
})

library(testthat)
library(data.table)

test_that(".check_numeric_columns() works as expected", {

    # Valid numeric columns
    dt1 <- data.table(a = 1:5, b = rnorm(5))
    expect_true(.check_numeric_columns(dt1, c("a", "b")))

    # Missing column
    dt2 <- data.table(x = 1:5)
    expect_error(
        .check_numeric_columns(dt2, c("x", "y")),
        "The following columns are not found in the data: y"
    )

    # Non-numeric column
    dt3 <- data.table(a = 1:5, b = letters[1:5])
    expect_error(
        .check_numeric_columns(dt3, c("a", "b")),
        "Non-numeric columns are in use.cols: b"
    )

    # Empty column list
    dt4 <- data.table(a = 1:5, b = 6:10)
    expect_true(.check_numeric_columns(dt4, character(0)))
})

test_that(".add_or_replace_column() adds or replaces columns with warning", {

    # Add a new column
    dt1 <- data.table(x = 1:3)
    result1 <- .add_or_replace_column(dt1, "y", dt1$x * 2)
    expect_true("y" %in% names(result1))
    expect_equal(result1$y, dt1$x * 2)

    # Replace existing column with warning
    dt2 <- data.table(x = 1:3, z = 3:1)
    expect_warning(
        result2 <- .add_or_replace_column(dt2, "z", dt2$x + 1),
        "Column 'z' already exists and will be replaced."
    )
    expect_equal(result2$z, dt2$x + 1)
})

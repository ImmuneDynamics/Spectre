test_that("make.colour.plot works for continuous colour", {
    dat <- data.table(
        X = rnorm(100),
        Y = rnorm(100),
        Value = rnorm(100)
    )
    p <- make.colour.plot(
        dat = dat,
        x.axis = "X",
        y.axis = "Y",
        col.axis = "Value",
        col.type = "continuous",
        save.to.disk = FALSE
    )
    expect_s3_class(p, "ggplot")
})

test_that("make.colour.plot works for factor colour with labels", {
    dat <- data.table(
        X = rnorm(100),
        Y = rnorm(100),
        Group = sample(letters[1:3], 100, replace = TRUE)
    )
    p <- make.colour.plot(
        dat = dat,
        x.axis = "X",
        y.axis = "Y",
        col.axis = "Group",
        col.type = "factor",
        add.label = TRUE,
        save.to.disk = FALSE
    )
    expect_s3_class(p, "ggplot")
})

test_that("make.colour.plot works for density plot (no col.axis)", {
    dat <- data.table(
        X = rnorm(100),
        Y = rnorm(100)
    )
    p <- make.colour.plot(
        dat = dat,
        x.axis = "X",
        y.axis = "Y",
        save.to.disk = FALSE
    )
    expect_s3_class(p, "ggplot")
})

test_that("make.colour.plot works with fast = TRUE", {
    dat <- data.table(
        X = rnorm(100),
        Y = rnorm(100),
        Group = sample(letters[1:3], 100, replace = TRUE)
    )
    p <- make.colour.plot(
        dat = dat,
        x.axis = "X",
        y.axis = "Y",
        col.axis = "Group",
        col.type = "factor",
        fast = TRUE,
        save.to.disk = FALSE
    )
    expect_s3_class(p, "ggplot")
})

test_that("make.colour.plot works with fast = TRUE for continuous colour", {
    dat <- data.table(
        X = rnorm(100),
        Y = rnorm(100),
        Value = rnorm(100)
    )
    p <- make.colour.plot(
        dat = dat,
        x.axis = "X",
        y.axis = "Y",
        col.axis = "Value",
        col.type = "continuous",
        fast = TRUE,
        save.to.disk = FALSE
    )
    expect_s3_class(p, "ggplot")
})

test_that("make.colour.plot works with different colour schemes", {
    dat <- data.table(
        X = rnorm(100),
        Y = rnorm(100),
        Value = rnorm(100),
        Group = sample(letters[1:3], 100, replace = TRUE)
    )
    # Test spectral (continuous)
    p1 <- make.colour.plot(
        dat = dat,
        x.axis = "X",
        y.axis = "Y",
        col.axis = "Value",
        col.type = "continuous",
        colours = "spectral",
        save.to.disk = FALSE
    )
    expect_s3_class(p1, "ggplot")

    # Test jet (continuous)
    p2 <- make.colour.plot(
        dat = dat,
        x.axis = "X",
        y.axis = "Y",
        col.axis = "Value",
        col.type = "continuous",
        colours = "jet",
        save.to.disk = FALSE
    )
    expect_s3_class(p2, "ggplot")

    # Test viridis (continuous)
    p3 <- make.colour.plot(
        dat = dat,
        x.axis = "X",
        y.axis = "Y",
        col.axis = "Value",
        col.type = "continuous",
        colours = "viridis",
        save.to.disk = FALSE
    )
    expect_s3_class(p3, "ggplot")

    # Test unlisted colour still work
    p4 <- make.colour.plot(
        dat = dat,
        x.axis = "X",
        y.axis = "Y",
        col.axis = "Group",
        col.type = "continuous",
        colours = "xxx",
        save.to.disk = FALSE
    )
    expect_s3_class(p4, "ggplot")
})
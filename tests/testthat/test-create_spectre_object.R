test_that("default create Spectre object works", {

    test_obj = create.spectre.object("cell_id", name = c("test"))
    
    expect_equal(test_obj@cell_id_col, "cell_id")
    expect_equal(test_obj@name, c("test"))
})

test_that("creating Spectre object without cell_id_col errors", {
    expect_error(
        create.spectre.object()
    )
})

test_that("creating Spectre object without name DO NOT errors", {
    expect_no_error(
        create.spectre.object("cell_id")
    )
})

# test adding new data

test_that("add new data works", {
    
    test_obj = create.spectre.object("cell_id")
    
    x = data.table(cell_id=seq(1,10), cd3=seq(1,10))
    test_obj = add.new.data(
        spectre_obj = test_obj, 
        dat = x, 
        dat_name = "cyto_test",
        metadata = list("description" = "test")
    )
    
    expect_equal(ncol(test_obj$cyto_test), 2)
    expect_equal(nrow(test_obj$cyto_test), 10)
    expect_equal(attributes(test_obj$cyto_test)$description, "test")
    
})

test_that("add new data with no metadata works", {

    test_obj = create.spectre.object("cell_id")

    x = data.table(cell_id=seq(1,10), cd3=seq(1,10))
    test_obj = add.new.data(
        spectre_obj = test_obj,
        dat = x,
        dat_name = "cyto_test"
    )

    expect_equal(ncol(test_obj$cyto_test), 2)
    expect_equal(nrow(test_obj$cyto_test), 10)
    expect_false("cyto_test" %in% names(attributes(test_obj$cyto_test)))


})

test_that("add new data works missing cell_id_col errors", {

    test_obj = create.spectre.object("cell_id")

    x = data.table(cd3=seq(1,10))
    expect_error(
        add.new.data(test_obj, x, "cyto_test")
    )
})

# test adding new metadata

test_that("add new metadata data works", {

    test_obj = create.spectre.object("cell_id")

    x = data.table(cell_id=seq(1,10), cd3=seq(1,10))
    test_obj = add.new.data(
        spectre_obj = test_obj,
        dat = x,
        dat_name = "cyto_test",
        metadata = list("description" = "test")
    )

    test_obj = add.new.metadata(
        spectre_obj = test_obj,
        metadata = list("test_descp" = c("test vector")),
        dataset_name = "cyto_test"
    )

    expect_true("description" %in% names(attributes(test_obj$cyto_test)))
    expect_true("test_descp" %in% names(attributes(test_obj$cyto_test)))
    expect_equal(attributes(test_obj$cyto_test)$description, "test")
    expect_equal(attributes(test_obj$cyto_test)$test_descp, c("test vector"))

})

test_that("add new metadata errors if the dataset doesn't exist", {

    test_obj = create.spectre.object("cell_id")

    x = data.table(cell_id=seq(1,10), cd3=seq(1,10))
    test_obj = add.new.data(
        spectre_obj = test_obj,
        dat = x,
        dat_name = "cyto_test",
        metadata = list("description" = "test")
    )

    expect_error(
        add.new.metadata(
            spectre_obj = test_obj,
            metadata = list("test_descp" = c("test vector")),
            dataset_name = "cyto_test2"
        )
    )

})



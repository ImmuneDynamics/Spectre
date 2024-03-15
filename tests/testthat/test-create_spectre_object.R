test_that("default create Spectre object works", {

    test_obj = create.spectre.object("cell_id")
    
    expect_equal(test_obj@cell_id_col, "cell_id")
})

test_that("no cell_id_col slot for creating Spectre object creation errors", {
    expect_error(
        create.spectre.object()
    )
})

test_that("add new data works", {
    
    test_obj = create.spectre.object("cell_id")
    
    x = data.table(cell_id=seq(1,10), cd3=seq(1,10))
    test_obj = add.new.data(test_obj, x, "cyto_test")
    
    expect_equal(ncol(test_obj$cyto_test), 2)
    expect_equal(nrow(test_obj$cyto_test), 10)
    
})

test_that("add new data works missing cell_id_col errors", {
    
    test_obj = create.spectre.object("cell_id")
    
    x = data.table(cd3=seq(1,10))
    expect_error(
        add.new.data(test_obj, x, "cyto_test")
    )
})


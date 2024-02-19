test_that("make autograph works", {
    dat <- data.table::data.table(
        Samples = c("Mock_01", "Mock_02", "Mock_03", "WNV_01", "WNV_02", "WNV_03"), 
        Group = c(rep("Mock",3), rep("WNV", 3)), 
        Tcells = c(20, 40, 30, 60, 70, 80), 
        Bcells = c(90, 95, 70, 20, 15, 30), 
        Batch = c(1, 2, 1, 3, 2, 1))
    
    expect_no_error(
        make.autograph(
            dat = dat, 
            x.axis = "Group",
            y.axis = "Tcells", 
            colour.by = "Batch", 
            colours = c("Black", "Red","Blue"), 
            y.axis.label = "Proportion",
            save_to_disk = FALSE)
    )
    
})

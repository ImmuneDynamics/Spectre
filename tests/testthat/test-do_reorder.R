test_that("do.reorder works", {
    dat <- Spectre::demo.clustered
    use.col <- 'FileName'
    
    new.order <- c('CNS_Mock_03.csv','CNS_Mock_04.csv','CNS_Mock_01.csv',
                   'CNS_Mock_02.csv','CNS_Mock_05.csv','CNS_Mock_06.csv',
                   'CNS_WNV_D7_01.csv','CNS_WNV_D7_02.csv','CNS_WNV_D7_03.csv',
                   'CNS_WNV_D7_04.csv','CNS_WNV_D7_05.csv','CNS_WNV_D7_06.csv')
    dat_reordered <- do.reorder(dat, use.col, new.order)
    
    actual_new_order <- unique(dat_reordered[[use.col]])
    
    for (i in seq(1, length(actual_new_order))) {
        expect_equal(
            as.character(actual_new_order[i]), 
            new.order[i]
        )
    }
})

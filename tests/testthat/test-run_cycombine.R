test_that("cycombine batch correction works", {
    dat = qs::qread(test_path("testdata", "artificial_batch_effect.qs"))
    
    dat = Spectre::do.subsample(dat, targets = rep(1000, 2), 
                                divide.by = "Batch")
    
    dat[, cell_id := paste0("Cell_", seq(nrow(dat)))]
    
    markers = c("CD45_chn", "CD48_chn", "CD117_chn", 
                "CD11b_chn", "SiglecF_chn", "NK11_chn", "B220_chn", 
                "CD8a_chn", "CD4_chn", "Ly6C_chn", "Ly6G_chn", "CD115_chn", 
                "CD3e_chn", "CD16.32_chn", "MHCII_chn")
    
    suppressWarnings(
        dat_corrected <- run.cycombine(
            dat = dat,
            cell_id_col = "cell_id",
            use_cols = markers,
            batch_col = "Batch",
            verbose = FALSE
        )
    )
    
    # just check there is a new element
    expect_true("cell_id" %in% names(dat_corrected))
    expect_true("Batch" %in% names(dat_corrected))
    expect_equal(nrow(dat), nrow(dat_corrected))
    
    # check the metadata
    expect_true("parameter" %in% names(attributes(dat_corrected)))
    expect_true("cycombine_extra_info" %in% names(attributes(dat_corrected)))
    
    # check the original data's column name has not been changed.
    expect_true("Batch" %in% names(dat_corrected))
    expect_false("batch" %in% names(dat_corrected))
})


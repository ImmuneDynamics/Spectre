test_that("rpca works", {
    dat = qs::qread(test_path("testdata", "artificial_batch_effect.qs"))
    
    # subsample because cca is so slow..
    dat = Spectre::do.subsample(dat, targets = rep(1000, 2), 
                                divide.by = "Batch")
    
    dat[, cell_id := paste0("Cell_", seq(nrow(dat)))]
    
    markers = c("CD45_chn", "CD48_chn", "CD117_chn", 
                "CD11b_chn", "SiglecF_chn", "NK11_chn", "B220_chn", 
                "CD8a_chn", "CD4_chn", "Ly6C_chn", "Ly6G_chn", "CD115_chn", 
                "CD3e_chn", "CD16.32_chn", "MHCII_chn")
    
    suppressWarnings(
        dat_corrected <- run.rpca(
            dat = dat,
            use_cols = markers,
            batch_col = "Batch",
            cell_id_col = "cell_id",
            verbose = FALSE,
            reference_batch = NULL
        )
    )
    
    # just check there is a new element
    expect_true(all(markers %in% names(dat_corrected)))
    expect_true("cell_id" %in% names(dat_corrected))
    expect_true("Batch" %in% names(dat_corrected))
    expect_equal(nrow(dat), nrow(dat_corrected))
})

test_that("rpca with reference works", {
    dat = qs::qread(test_path("testdata", "artificial_batch_effect.qs"))
    
    # subsample because cca is so slow..
    dat = Spectre::do.subsample(dat, targets = rep(1000, 2), 
                                divide.by = "Batch")
    
    dat[, cell_id := paste0("Cell_", seq(nrow(dat)))]
    
    markers = c("CD45_chn", "CD48_chn", "CD117_chn", 
                "CD11b_chn", "SiglecF_chn", "NK11_chn", "B220_chn", 
                "CD8a_chn", "CD4_chn", "Ly6C_chn", "Ly6G_chn", "CD115_chn", 
                "CD3e_chn", "CD16.32_chn", "MHCII_chn")
    
    suppressWarnings(
        dat_corrected <- run.rpca(
            dat = dat,
            use_cols = markers,
            batch_col = "Batch",
            cell_id_col = "cell_id",
            verbose = FALSE,
            reference_batch = "A"
        )
    )
    
    # just check there is a new element
    expect_true(all(markers %in% names(dat_corrected)))
    expect_true("cell_id" %in% names(dat_corrected))
    expect_true("Batch" %in% names(dat_corrected))
    expect_equal(nrow(dat), nrow(dat_corrected))
})



test_that("harmony with pca works", {
    dat = qs::qread(test_path("testdata", "artificial_batch_effect.qs"))
    
    dat[, cell_id := paste0("Cell_", seq(nrow(dat)))]
    
    markers = c("CD45_chn", "CD48_chn", "CD117_chn", 
                "CD11b_chn", "SiglecF_chn", "NK11_chn", "B220_chn", 
                "CD8a_chn", "CD4_chn", "Ly6C_chn", "Ly6G_chn", "CD115_chn", 
                "CD3e_chn", "CD16.32_chn", "MHCII_chn")
    
    suppressWarnings(
        res <- run.harmony(
            dat = dat,
            cell_id_col = "cell_id",
            use_cols = markers,
            batch_col = "Batch",
            verbose = FALSE
        )
    )
    
    # just check there is a new element
    expect_true("cell_id" %in% names(res))
    expect_true("Batch" %in% names(res))
    expect_equal(nrow(dat), nrow(res))
    expect_true(all(paste0("PC_", seq(length(markers) - 1)) %in% names(res)))
    
    # check the metadata
    expect_true("harmony_parameter" %in% names(attributes(res)))
})

test_that("harmony NO pca works", {
    dat = qs::qread(test_path("testdata", "artificial_batch_effect.qs"))
    
    dat[, cell_id := paste0("Cell_", seq(nrow(dat)))]
    
    markers = c("CD45_chn", "CD48_chn", "CD117_chn", 
                "CD11b_chn", "SiglecF_chn", "NK11_chn", "B220_chn", 
                "CD8a_chn", "CD4_chn", "Ly6C_chn", "Ly6G_chn", "CD115_chn", 
                "CD3e_chn", "CD16.32_chn", "MHCII_chn")
    
    suppressWarnings(
        res <- run.harmony(
            dat = dat,
            cell_id_col = "cell_id",
            use_cols = markers,
            batch_col = "Batch",
            verbose = FALSE,
            do_pca = FALSE
        )
    )
    
    # just check there is a new element
    expect_true("cell_id" %in% names(res))
    expect_true("Batch" %in% names(res))
    expect_equal(nrow(dat), nrow(res))
    expect_true(all(markers %in% names(res)))
    
    # check the metadata
    expect_true("harmony_parameter" %in% names(attributes(res)))
})

test_that("return_object works", {
    dat = qs::qread(test_path("testdata", "artificial_batch_effect.qs"))
    
    dat[, cell_id := paste0("Cell_", seq(nrow(dat)))]
    
    markers = c("CD45_chn", "CD48_chn", "CD117_chn", 
                "CD11b_chn", "SiglecF_chn", "NK11_chn", "B220_chn", 
                "CD8a_chn", "CD4_chn", "Ly6C_chn", "Ly6G_chn", "CD115_chn", 
                "CD3e_chn", "CD16.32_chn", "MHCII_chn")
    
    suppressWarnings(
        res <- run.harmony(
            dat = dat,
            cell_id_col = "cell_id",
            use_cols = markers,
            batch_col = "Batch",
            verbose = FALSE,
            return_object = TRUE
        )
    )
    
    expect_true("harmony_object" %in% names(attributes(res)))
})




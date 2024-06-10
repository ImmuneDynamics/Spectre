test_that("rpca batch correction works", {
    dat_raw = Spectre::demo.batches
    
    # subsample for speedy execution
    dat_raw = Spectre::do.subsample(dat_raw, targets = rep(1000, 2), 
                                    divide.by = "Batch")
    
    dat_raw[, cell_id := paste0("Cell_", seq(nrow(dat_raw)))]
    dat = create.spectre.object(cell_id_col = "cell_id")
    dat = add.new.data(spectre_obj = dat, dat = dat_raw, "cyto_batch")
    
    markers = c("CD45_chn", "CD48_chn", "CD117_chn", 
                "CD11b_chn", "SiglecF_chn", "NK11_chn", "B220_chn", 
                "CD8a_chn", "CD4_chn", "Ly6C_chn", "Ly6G_chn", "CD115_chn", 
                "CD3e_chn", "CD16.32_chn", "MHCII_chn")
    
    suppressWarnings(
        dat <- run.rpca.batch.correction(
            dat = dat,
            data_source = "cyto_batch",
            output_name = "cyto_batch_corrected",
            use_cols = markers,
            batch_col = "Batch",
            verbose = FALSE,
            reference_batch = NULL
        )
    )
    
    # just check there is a new element
    expect_true(is.data.table(dat$cyto_batch_corrected))
    expect_true(all(markers %in% names(dat$cyto_batch_corrected)))
    expect_true("cell_id" %in% names(dat$cyto_batch_corrected))
    expect_true("Batch" %in% names(dat$cyto_batch_corrected))
    expect_equal(nrow(dat$cyto_batch), nrow(dat$cyto_batch_corrected))
})

test_that("rpca batch correction with reference works", {
    dat_raw = Spectre::demo.batches
    
    # subsample for speedy execution
    dat_raw = Spectre::do.subsample(dat_raw, targets = rep(1000, 2), 
                                    divide.by = "Batch")
    
    
    dat_raw[, cell_id := paste0("Cell_", seq(nrow(dat_raw)))]
    dat = create.spectre.object(cell_id_col = "cell_id")
    dat = add.new.data(spectre_obj = dat, dat = dat_raw, "cyto_batch")
    
    markers = c("CD45_chn", "CD48_chn", "CD117_chn", 
                "CD11b_chn", "SiglecF_chn", "NK11_chn", "B220_chn", 
                "CD8a_chn", "CD4_chn", "Ly6C_chn", "Ly6G_chn", "CD115_chn", 
                "CD3e_chn", "CD16.32_chn", "MHCII_chn")
    
    suppressWarnings(
        dat <- run.rpca.batch.correction(
            dat = dat,
            data_source = "cyto_batch",
            output_name = "cyto_batch_corrected",
            use_cols = markers,
            batch_col = "Batch",
            verbose = FALSE,
            reference_batch = "A"
        )
    )
    
    # just check there is a new element
    expect_true(is.data.table(dat$cyto_batch_corrected))
    expect_true(all(markers %in% names(dat$cyto_batch_corrected)))
    expect_true("cell_id" %in% names(dat$cyto_batch_corrected))
    expect_true("Batch" %in% names(dat$cyto_batch_corrected))
    expect_equal(nrow(dat$cyto_batch), nrow(dat$cyto_batch_corrected))
})


test_that("rpca batch correction for non-spectre object fails", {
    dat_raw = Spectre::demo.batches
    
    # subsample for speedy execution
    dat_raw = Spectre::do.subsample(dat_raw, targets = rep(1000, 2), 
                                    divide.by = "Batch")
    
    
    markers = c("CD45_chn", "CD48_chn", "CD117_chn", 
                "CD11b_chn", "SiglecF_chn", "NK11_chn", "B220_chn", 
                "CD8a_chn", "CD4_chn", "Ly6C_chn", "Ly6G_chn", "CD115_chn", 
                "CD3e_chn", "CD16.32_chn", "MHCII_chn")
    
    expect_error(
        suppressWarnings(
            run.rpca.batch.correction(
                dat = dat_raw,
                data_source = "cyto_batch",
                output_name = "cyto_batch_corrected",
                use_cols = markers,
                batch_col = "Batch"
            )
        )
    )
})


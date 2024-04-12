test_that("rpca label_transfer infer common markers works", {
    dat_raw = Spectre::demo.batches
    
    # subsample for speedy execution
    dat_raw = Spectre::do.subsample(dat_raw, targets = rep(1000, 2), 
                                    divide.by = "Batch")
    
    dat_raw[, cell_id := paste0("Cell_", seq(nrow(dat_raw)))]
    
    dat_raw <- split(dat_raw, dat_raw$Batch)
    
    # get rid of all the non markers columns
    markers = c("CD45_chn", "CD48_chn", "CD117_chn",
                "CD11b_chn", "SiglecF_chn", "NK11_chn", "B220_chn",
                "CD8a_chn", "CD4_chn", "Ly6C_chn", "Ly6G_chn", "CD115_chn",
                "CD3e_chn", "CD16.32_chn", "MHCII_chn", "Population", "cell_id")
    
    dat_raw <- lapply(dat_raw, function(d) d[, markers, with = FALSE])
    dat_raw$A[, Population := NULL]
    
    dat = create.spectre.object(cell_id_col = "cell_id")
    dat = add.new.data(spectre_obj = dat, dat = dat_raw$A, "cyto")
    dat = add.new.data(spectre_obj = dat, dat = dat_raw$B, "cite")
    
    suppressWarnings(
        dat <- run.rpca.label.transfer(
            dat = dat,
            cytometry_data_source = "cyto",
            citeseq_data_source = "cite",
            output_name = "labels",
            celltype_label_col = "Population",
            verbose = FALSE
        )
    )
    
    # just check there is a new element
    expect_true("predicted_cell_type" %in% names(dat$labels))
    expect_true("cell_id" %in% names(dat$labels))
    expect_true("rpca_prediction_score" %in% names(attributes(dat$labels)))
    expect_equal(nrow(dat$labels), nrow(dat$cyto))
})

test_that("rpca label_transfer valid use_cols works", {
    dat_raw = Spectre::demo.batches
    
    # subsample for speedy execution
    dat_raw = Spectre::do.subsample(dat_raw, targets = rep(1000, 2), 
                                    divide.by = "Batch")
    
    dat_raw[, cell_id := paste0("Cell_", seq(nrow(dat_raw)))]
    
    dat_raw <- split(dat_raw, dat_raw$Batch)
    
    # get rid of all the non markers columns
    markers = c("CD45_chn", "CD48_chn", "CD117_chn",
                "CD11b_chn", "SiglecF_chn", "NK11_chn", "B220_chn",
                "CD8a_chn", "CD4_chn", "Ly6C_chn", "Ly6G_chn", "CD115_chn",
                "CD3e_chn", "CD16.32_chn", "MHCII_chn", "Population", "cell_id")
    
    dat_raw <- lapply(dat_raw, function(d) d[, markers, with = FALSE])
    dat_raw$A[, Population := NULL]
    
    dat = create.spectre.object(cell_id_col = "cell_id")
    dat = add.new.data(spectre_obj = dat, dat = dat_raw$A, "cyto")
    dat = add.new.data(spectre_obj = dat, dat = dat_raw$B, "cite")
    
    suppressWarnings(
        dat <- run.rpca.label.transfer(
            dat = dat,
            cytometry_data_source = "cyto",
            citeseq_data_source = "cite",
            output_name = "labels",
            celltype_label_col = "Population",
            verbose = FALSE,
            use_cols = c("CD45_chn", "CD48_chn", "CD117_chn",
                         "CD11b_chn")
        )
    )
    
    # just check there is a new element
    expect_true("predicted_cell_type" %in% names(dat$labels))
    expect_true("cell_id" %in% names(dat$labels))
    expect_true("rpca_prediction_score" %in% names(attributes(dat$labels)))
    expect_equal(nrow(dat$labels), nrow(dat$cyto))
})

test_that("rpca label_transfer with semi valid use_cols works", {
    dat_raw = Spectre::demo.batches
    
    # subsample for speedy execution
    dat_raw = Spectre::do.subsample(dat_raw, targets = rep(1000, 2), 
                                    divide.by = "Batch")
    
    dat_raw[, cell_id := paste0("Cell_", seq(nrow(dat_raw)))]
    
    dat_raw <- split(dat_raw, dat_raw$Batch)
    
    # get rid of all the non markers columns
    markers = c("CD45_chn", "CD48_chn", "CD117_chn",
                "CD11b_chn", "SiglecF_chn", "NK11_chn", "B220_chn",
                "CD8a_chn", "CD4_chn", "Ly6C_chn", "Ly6G_chn", "CD115_chn",
                "CD3e_chn", "CD16.32_chn", "MHCII_chn", "Population", "cell_id")
    
    dat_raw <- lapply(dat_raw, function(d) d[, markers, with = FALSE])
    dat_raw$A[, Population := NULL]
    
    dat = create.spectre.object(cell_id_col = "cell_id")
    dat = add.new.data(spectre_obj = dat, dat = dat_raw$A, "cyto")
    dat = add.new.data(spectre_obj = dat, dat = dat_raw$B, "cite")
    
    suppressWarnings(
        dat <- run.rpca.label.transfer(
            dat = dat,
            cytometry_data_source = "cyto",
            citeseq_data_source = "cite",
            output_name = "labels",
            celltype_label_col = "Population",
            verbose = FALSE,
            use_cols = c("CD45_chn", "CD48_chn", "CD117_chn",
                         "X")
        )
    )
    
    # just check there is a new element
    expect_true("predicted_cell_type" %in% names(dat$labels))
    expect_true("cell_id" %in% names(dat$labels))
    expect_true("rpca_prediction_score" %in% names(attributes(dat$labels)))
    expect_equal(nrow(dat$labels), nrow(dat$cyto))
})


test_that("rpca no common markers fails", {
    dat_cite <- data.table(
        cell_id = paste0("cell_", seq(1, 10)),
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2),
        labels = paste0("pop_", seq(1, 10))
    )
    
    dat_cyto <- data.table(
        cell_id = paste0("cell_", seq(1, 10)),
        marker3 = rnorm(10, 3)
    )
    
    dat = create.spectre.object(cell_id_col = "cell_id")
    dat = add.new.data(spectre_obj = dat, dat = dat_cyto, "cyto")
    dat = add.new.data(spectre_obj = dat, dat = dat_cite, "cite")
    
    expect_error(
        suppressWarnings(
            run.rpca.label.transfer(
                dat = dat,
                cytometry_data_source = "cyto",
                citeseq_data_source = "cite",
                output_name = "labels",
                celltype_label_col = "Population",
                verbose = FALSE
            )
        )
    )
    
    expect_error(
        suppressWarnings(
            run.rpca.label.transfer(
                dat = dat,
                cytometry_data_source = "cyto",
                citeseq_data_source = "cite",
                output_name = "labels",
                celltype_label_col = "labels",
                verbose = FALSE,
                use_cols = c("marker_3")
            )
        )
    )
    
    expect_error(
        suppressWarnings(
            run.rpca.label.transfer(
                dat = dat,
                cytometry_data_source = "cyto",
                citeseq_data_source = "cite",
                output_name = "labels",
                celltype_label_col = "labels",
                verbose = FALSE,
                use_cols = c("marker3")
            )
        )
    )
    
    expect_error(
        suppressWarnings(
            run.rpca.label.transfer(
                dat = dat,
                cytometry_data_source = "cyto",
                citeseq_data_source = "cite",
                output_name = "labels",
                celltype_label_col = "labels",
                verbose = FALSE,
                use_cols = c("marker1")
            )
        )
    )
})


test_that("CSV files are read into a list", {
    dat <- read.files(
        file.loc = test_path("testdata", "csv_fcs_files"),
        verbose = FALSE
    )
    expect_equal(typeof(dat), "list")
    expect_equal(length(dat), 2)
    
    file_names <- c("test_Mock_01_subsampled_1", "test_Mock_01_subsampled_2")
    
    expect_equal(names(dat), c(file_names))

    # Test content
    
    for (i in seq(1, length(file_names))) {
        f <- file_names[[i]]
        expect_equal(nrow(dat[[f]]), 500)
        expect_equal(ncol(dat[[f]]), 10)
        expect_equal(
            names(dat[[f]]),
            c(
                "PECy5-5 CD3e", "PECy7 CD16-32", "DL800 Ly6G",
                "AF700 CD45", "APCCy7 CD48", "BUV395 CD11b",
                "BUV737 B220", "BV605 Ly6C", "FileName", "FileNo"
            )
        )
        expect_equal(
            dat[[f]]$FileName,
            rep(f, 500)
        )
        
        expect_equal(
            dat[[f]]$FileNo,
            rep(i, 500)
        )
    }
})

test_that("CSV files are read into a list", {
    dat <- read.files(
        file.loc = test_path("testdata", "csv_fcs_files"),
        verbose = FALSE
    )
    expect_equal(typeof(dat), "list")
    expect_equal(length(dat), 2)
    
    file_names <- c("test_Mock_01_subsampled_1", "test_Mock_01_subsampled_2")
    
    expect_equal(names(dat), c(file_names))

    # Test content
    
    for (i in seq(1, length(file_names))) {
        f <- file_names[[i]]
        expect_equal(nrow(dat[[f]]), 500)
        expect_equal(ncol(dat[[f]]), 10)
        expect_equal(
            names(dat[[f]]),
            c(
                "PECy5-5 CD3e", "PECy7 CD16-32", "DL800 Ly6G",
                "AF700 CD45", "APCCy7 CD48", "BUV395 CD11b",
                "BUV737 B220", "BV605 Ly6C", "FileName", "FileNo"
            )
        )
        expect_equal(
            dat[[f]]$FileName,
            rep(f, 500)
        )
        
        expect_equal(
            dat[[f]]$FileNo,
            rep(i, 500)
        )
    }
})

test_that("CSV files are read into a SpectreObject", {
    dat <- read.files(
        file.loc = test_path("testdata", "csv_fcs_files"),
        as_spectre_object = TRUE,
        verbose = FALSE
    )
    
    expect_equal(nrow(dat@cytometry_data), 1000)
    expect_equal(ncol(dat@cytometry_data), 10)
    
    expect_equal(
        names(dat@cytometry_data),
        c(
            "PECy5-5 CD3e", "PECy7 CD16-32", "DL800 Ly6G",
            "AF700 CD45", "APCCy7 CD48", "BUV395 CD11b",
            "BUV737 B220", "BV605 Ly6C", "FileName", "FileNo"
        )
    )
    
    expect_equal(
        dat@cytometry_data$FileName,
        c(
            rep("test_Mock_01_subsampled_1", 500),
            rep("test_Mock_01_subsampled_2", 500)
        )
    )
    
    expect_equal(
        dat@cytometry_data$FileNo,
        c(
            rep(1, 500),
            rep(2, 500)
        )
    )
    
})

test_that("FCS files are read into SpectreObject", {
    dat <- read.files(
        file.loc = test_path("testdata", "csv_fcs_files"),
        file.type = ".fcs",
        as_spectre_object = TRUE,
        verbose = FALSE
    )
    

    # Test content
    expect_equal(nrow(dat@cytometry_data), 1000)
    expect_equal(ncol(dat@cytometry_data), 10)

    # Weird naming of markers if read FCS..why bother?
    
    expect_equal(
        names(dat@cytometry_data), 
        c(
            "PECy5-5 CD3e_column PECy5-5 CD3e from dataset",
            "PECy7 CD16-32_column PECy7 CD16-32 from dataset",
            "DL800 Ly6G_column DL800 Ly6G from dataset",
            "AF700 CD45_column AF700 CD45 from dataset",
            "APCCy7 CD48_column APCCy7 CD48 from dataset",
            "BUV395 CD11b_column BUV395 CD11b from dataset",
            "BUV737 B220_column BUV737 B220 from dataset",
            "BV605 Ly6C_column BV605 Ly6C from dataset",
            "FileName",
            "FileNo"
        )    
    )
})

test_that("FCS files are read into list", {
    dat <- read.files(
        file.loc = test_path("testdata", "csv_fcs_files"),
        file.type = ".fcs",
        verbose = FALSE
    )
    expect_equal(typeof(dat), "list")
    expect_equal(length(dat), 1)
    expect_equal(names(dat), c("test_Mock_01_subsampled"))
    
    # Test content
    dat <- dat[["test_Mock_01_subsampled"]]
    expect_equal(nrow(dat), 1000)
    expect_equal(ncol(dat), 10)
    
    # Weird naming of markers if read FCS..why bother?
    expected_headers <- c(
        "PECy5-5 CD3e_column PECy5-5 CD3e from dataset",
        "PECy7 CD16-32_column PECy7 CD16-32 from dataset",
        "DL800 Ly6G_column DL800 Ly6G from dataset",
        "AF700 CD45_column AF700 CD45 from dataset",
        "APCCy7 CD48_column APCCy7 CD48 from dataset",
        "BUV395 CD11b_column BUV395 CD11b from dataset",
        "BUV737 B220_column BUV737 B220 from dataset",
        "BV605 Ly6C_column BV605 Ly6C from dataset",
        "FileName",
        "FileNo"
    )
    
    expect_equal(names(dat), expected_headers)
})

test_that("read files without embed filenames", {
    dat_csv <- read.files(
        file.loc = test_path("testdata", "csv_fcs_files"),
        do.embed.file.names = FALSE,
        verbose = FALSE
    )
    dat_fcs <- read.files(
        file.loc = test_path("testdata", "csv_fcs_files"),
        do.embed.file.names = FALSE,
        file.type = '.fcs',
        verbose = FALSE
    )
    dat_obj_csv <- read.files(
        file.loc = test_path("testdata", "csv_fcs_files"),
        do.embed.file.names = FALSE,
        as_spectre_object = TRUE,
        verbose = FALSE
    )
    dat_obj_fcs <- read.files(
        file.loc = test_path("testdata", "csv_fcs_files"),
        do.embed.file.names = FALSE,
        as_spectre_object = TRUE,
        file.type = ".fcs",
        verbose = FALSE
    )
    
    file_names <- c("test_Mock_01_subsampled_1", "test_Mock_01_subsampled_2")

    # Test content csv
    for (i in seq(1, length(file_names))) {
        f <- file_names[[i]]
        
        expect_equal(
            names(dat_csv[[f]]),
            c(
                "PECy5-5 CD3e", "PECy7 CD16-32", "DL800 Ly6G",
                "AF700 CD45", "APCCy7 CD48", "BUV395 CD11b",
                "BUV737 B220", "BV605 Ly6C"
            )
        )
    }
    
    # Test content fcs
    expect_equal(
        names(dat_fcs[["test_Mock_01_subsampled"]]),
        c(
            "PECy5-5 CD3e_column PECy5-5 CD3e from dataset",
            "PECy7 CD16-32_column PECy7 CD16-32 from dataset",
            "DL800 Ly6G_column DL800 Ly6G from dataset",
            "AF700 CD45_column AF700 CD45 from dataset",
            "APCCy7 CD48_column APCCy7 CD48 from dataset",
            "BUV395 CD11b_column BUV395 CD11b from dataset",
            "BUV737 B220_column BUV737 B220 from dataset",
            "BV605 Ly6C_column BV605 Ly6C from dataset"
        )
    )
    
    # Test content SpectreObject
    expect_equal(
        names(dat_obj_csv@cytometry_data),
        c(
            "PECy5-5 CD3e", "PECy7 CD16-32", "DL800 Ly6G",
            "AF700 CD45", "APCCy7 CD48", "BUV395 CD11b",
            "BUV737 B220", "BV605 Ly6C"
        )
    )
    
    expect_equal(
        names(dat_obj_fcs@cytometry_data),
        c(
            "PECy5-5 CD3e_column PECy5-5 CD3e from dataset",
            "PECy7 CD16-32_column PECy7 CD16-32 from dataset",
            "DL800 Ly6G_column DL800 Ly6G from dataset",
            "AF700 CD45_column AF700 CD45 from dataset",
            "APCCy7 CD48_column APCCy7 CD48 from dataset",
            "BUV395 CD11b_column BUV395 CD11b from dataset",
            "BUV737 B220_column BUV737 B220 from dataset",
            "BV605 Ly6C_column BV605 Ly6C from dataset"
        )
    )
})

test_that("missing file.loc result in error", {
    expect_error(read.files())
})

test_that("missing files result in error", {
    expect_error(read.files(file.loc = test_path("testdata")))
    
    expect_error(
        read.files(
            file.loc = test_path("testdata"),
            file.type = ".fcs",
            verbose = FALSE
        )
    )
})

test_that("unsupported file type result in error", {
    expect_error(read.files(".txt"))
})





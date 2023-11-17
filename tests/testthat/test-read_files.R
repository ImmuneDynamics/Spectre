# TODO: need to be updated for v2
#' #' Unit test for reading files. Data are stored in inst/testdata
#' 
#' test_that("CSV files are read", {
#'     dir_loc <- system.file("testdata", package = "Spectre")
#'     dat <- read.files.to.list(file.loc = dir_loc)
#'     expect_equal(typeof(dat), "list")
#'     expect_equal(length(dat), 1)
#'     expect_equal(names(dat), c("test_Mock_01_subsampled"))
#' 
#'     # Test content
#'     dat <- dat[["test_Mock_01_subsampled"]]
#'     expect_equal(nrow(dat), 1000)
#'     expect_equal(ncol(dat), 10)
#'     expect_equal(
#'         names(dat),
#'         c(
#'             "PECy5-5 CD3e", "PECy7 CD16-32", "DL800 Ly6G",
#'             "AF700 CD45", "APCCy7 CD48", "BUV395 CD11b",
#'             "BUV737 B220", "BV605 Ly6C", "FileName", "FileNo"
#'         )
#'     )
#' })
#' 
#' test_that("FCS files are read", {
#'     dir_loc <- system.file("testdata", package = "Spectre")
#'     dat <- read.files.to.list(file.loc = dir_loc, file.type = ".fcs")
#'     expect_equal(typeof(dat), "list")
#'     expect_equal(length(dat), 1)
#'     expect_equal(names(dat), c("test_Mock_01_subsampled"))
#' 
#'     # Test content
#'     dat <- dat[["test_Mock_01_subsampled"]]
#'     expect_equal(nrow(dat), 1000)
#'     expect_equal(ncol(dat), 10)
#' 
#'     # Weird naming of markers if read FCS..why bother?
#'     expected_headers <- c(
#'         "PECy5-5 CD3e_column PECy5-5 CD3e from dataset",
#'         "PECy7 CD16-32_column PECy7 CD16-32 from dataset",
#'         "DL800 Ly6G_column DL800 Ly6G from dataset",
#'         "AF700 CD45_column AF700 CD45 from dataset",
#'         "APCCy7 CD48_column APCCy7 CD48 from dataset",
#'         "BUV395 CD11b_column BUV395 CD11b from dataset",
#'         "BUV737 B220_column BUV737 B220 from dataset",
#'         "BV605 Ly6C_column BV605 Ly6C from dataset",
#'         "FileName",
#'         "FileNo"
#'     )
#' 
#'     expect_equal(names(dat), expected_headers)
#' })
#' 
#' test_that("FCS files with truncate set to false are read", {
#'     dir_loc <- system.file("testdata", package = "Spectre")
#'     dat <- read.files.to.list(file.loc = dir_loc, file.type = ".fcs", truncate_max_range = FALSE)
#'     expect_equal(typeof(dat), "list")
#'     expect_equal(length(dat), 1)
#'     expect_equal(names(dat), c("test_Mock_01_subsampled"))
#' 
#'     # Test content
#'     dat <- dat[["test_Mock_01_subsampled"]]
#'     expect_equal(nrow(dat), 1000)
#'     expect_equal(ncol(dat), 10)
#' 
#'     # Weird naming of markers if read FCS..why bother?
#'     expected_headers <- c(
#'         "PECy5-5 CD3e_column PECy5-5 CD3e from dataset",
#'         "PECy7 CD16-32_column PECy7 CD16-32 from dataset",
#'         "DL800 Ly6G_column DL800 Ly6G from dataset",
#'         "AF700 CD45_column AF700 CD45 from dataset",
#'         "APCCy7 CD48_column APCCy7 CD48 from dataset",
#'         "BUV395 CD11b_column BUV395 CD11b from dataset",
#'         "BUV737 B220_column BUV737 B220 from dataset",
#'         "BV605 Ly6C_column BV605 Ly6C from dataset",
#'         "FileName",
#'         "FileNo"
#'     )
#' 
#'     expect_equal(names(dat), expected_headers)
#' })

#' run.ruv - Run the the CytofRUV alignment method.
#'
#' This function allows you to run the CytofRUV alignment method.
#'
#' @usage run.align()
#'
#' @param dat NO DEFAULT.A data.table consisting of the data (both reference and taret) you wish to align
#' @param sample.col NO DEFAULT. Column name of the data.table that contains sample names
#' @param batch.col NO DEFAULT. Column name of the data.table that contains batch labels
#' @param align.cols NO DEFAULT. Vector of column names to align.
#' @param ref.samples NO DEFAULT. Character vector of sample names that are used as reference data.
#' @param cluster.col DEFAULT = NULL. Column name of the data.table that contains clusters. If NULL, will align data without clusters.
#' @param k DEFAULT = K. K value for alignment
#' @param dir DEFAULT = getwd(). Directory where the work is being done
#' @param append.name DEFAULT = '_ruv'. Text to be appended to end of each new column containing aligned data
#'
#' @return Returns a data.table with aligned data added in new columns.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
#'
#'
#' @import data.table
#'

run.ruv <- function(dat,
                    sample.col,
                    batch.col,
                    align.cols,
                    ref.samples,
                    cluster.col = NULL,
                    k = 5,
                    dir = getwd(),
                    append.name = "_ruv") {

  ### Packages

  if (!is.element("data.table", installed.packages()[, 1])) stop("data.table is required but not installed")
  if (!is.element("CytofRUV", installed.packages()[, 1])) stop("CytofRUV is required but not installed")

  require(data.table)
  require(CytofRUV)

  ### Test data

  # dat <- target.dat
  #
  # lst <- c("Mock01_A", "Mock02_A", "Mock03_A",
  #          "Mock04_B", "Mock05_B", "Mock06_B",
  #
  #          "WNV01_A", "WNV02_A", "WNV03_A",
  #          "WNV04_B", "WNV05_B", "WNV06_B")
  #
  # tb <- data.table("FileName" = sort(unique(dat$FileName)),
  #                  "NewName" = lst)
  #
  # dat <- do.add.cols(dat, "FileName", tb, "FileName")
  # setorderv(dat, 'NewName')
  # setorderv(dat, 'Batches')
  # dat
  #
  # # dat$cluster <- 1
  #
  # sample.col <- "Sample"
  # batch.col <- "Batches"
  #
  # align.cols
  # ref.samples <- c("CNS_WNV_D7_01", "CNS_WNV_D7_04")
  #
  # #cluster.col = 'cluster' # if NULL
  # cluster.col <- NULL
  #
  # k = 5
  # seed = 1234
  #
  # dir = getwd()
  # append.name = '_ruv'
  #
  # dat

  ### Prepare functions

  run_RUVIII <- function(data, norm_clusters, k, rep_samples) {
    raw_Y <- as.matrix(data[3:ncol(data)])
    # Standardise the input and then compensate output
    col_means <- colMeans(raw_Y)
    col_sds <- apply(raw_Y, 2, function(x) sd(x))

    for (i in 1:ncol(raw_Y)) {
      raw_Y[, i] <- (raw_Y[, i] - col_means[i]) / col_sds[i]
    }
    # Run the actual RUVIII
    res_mat <- make_residual_mat_several_rep(raw_Y, data$cluster, norm_clusters, data$sample, rep_samples)
    norm_Y <- fastRUVIII(Y = raw_Y, M, ctl = c(1:ncol(raw_Y)), res_mat = res_mat, k = k)$newY

    for (i in 1:ncol(norm_Y)) {
      norm_Y[, i] <- norm_Y[, i] * col_sds[i] + col_means[i]
    }

    return(norm_Y)
  }

  make_residual_mat_several_rep <- function(data, clusters, norm_clus_list, samples, rep_samples_list) {
    # make_residual_mat(raw_Y,data$cluster, norm_clusters,norm_clusters_second,data$sample,rep_samples,second_rep_samples)
    res_mat <- data
    res_mat[] <- 0
    for (r in 1:length(rep_samples_list)) {
      norm_clus <- norm_clus_list[[r]]
      rep_samples <- rep_samples_list[[r]]
      mean_pseud <- matrix(nrow = length(norm_clus), ncol = dim(data)[2])
      for (i in 1:length(norm_clus)) {
        tmp <- ((clusters == norm_clus[i]) & (samples %in% rep_samples))
        mean_pseud[i, ] <- colMeans(data[tmp, ])
        res_mat[tmp, ] <- t(apply(data[tmp, ], 1, function(x) x - mean_pseud[i, ]))
      }
    }
    return(res_mat)
  }

  fastRUVIII <- function(Y, M, ctl, res_mat, k = NULL, eta = NULL, average = FALSE, fullalpha = NULL) {
    # Assumes good input
    if (!(k > 0)) stop("Bad input - read the documentation")
    Y <- ruv::RUV1(Y, eta, ctl)
    m <- nrow(Y)
    # Y0 = fast_residop(Y, M)
    Y0 <- res_mat
    fullalpha <- diag(rsvd::rsvd(Y0)$d) %*% t(rsvd::rsvd(Y0)$v)
    alpha <- fullalpha[1:k, , drop = FALSE]
    ac <- alpha[, ctl, drop = FALSE]
    W <- Y[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
    newY <- Y - W %*% alpha
    return(list(newY = newY, fullalpha = fullalpha))
  }

  fast_residop <- function(A, B) {
    return(A - B %*% solve(t(B) %*% B) %*% (t(B) %*% A))
  }

  ### Checks

  ## Clusters present in all batches

  if (is.null(cluster.col)) {
    dat$`RUV_CLUST` <- 1
    cluster.col <- "RUV_CLUST"
  }

  ## Clusters present? If not, all = 1

  all.clusters <- unique(dat[[cluster.col]])

  for (i in unique(dat[[batch.col]])) {
    # i <- unique(dat[[batch.col]])[[1]]
    temp <- dat[dat[[batch.col]] == i, cluster.col, with = FALSE]
    temp <- temp[[1]]
    temp <- unique(temp)
    temp

    if (sort(temp) != sort(all.clusters)) {
      stop("Some clusters are missing from certain batches. For RUV to work, clusters need to be represented in each batch")
    }
  }

  ### Cluster data preparation

  start.dat <- dat

  clusters <- sort(unique(dat[[cluster.col]]))
  RUV_CLUSTER_NUMBERS <- c(1:length(clusters))

  clusters.tb <- data.table(
    "clusters" = clusters,
    "RUV_CLUSTER_NUMBERS" = RUV_CLUSTER_NUMBERS
  )
  dat <- do.add.cols(dat, cluster.col, clusters.tb, "clusters")

  ### Batch data preparation

  batch.names <- sort(unique(dat[[batch.col]]))
  RUV_BATCH_NUMBERS <- c(1:length(batch.names))

  batch.tb <- data.table(
    "batch.names" = batch.names,
    "RUV_BATCH_NUMBERS" = RUV_BATCH_NUMBERS
  )

  dat <- do.add.cols(dat, batch.col, batch.tb, "batch.names")

  ### Sample data preparation

  sample.names <- sort(unique(dat[[sample.col]]))
  sample.numbers <- paste0("Sample", c(1:length(sample.names)))

  sample.tb <- data.table(
    "sample.names" = sample.names,
    "RUV_SAMPLES" = sample.numbers
  )


  dat <- do.add.cols(dat, sample.col, sample.tb, "sample.names")
  dat$RUV_SAMPLES <- paste0(dat$RUV_SAMPLES, "_", dat$RUV_BATCH_NUMBERS)

  ### Sample batch matching

  sample.batches <- list()
  for (i in c(1:length(sample.names))) {
    # i <- 1
    a <- sample.names[[i]]
    temp <- dat[dat[[sample.col]] == a, "RUV_BATCH_NUMBERS", with = FALSE]
    temp <- temp[1, ]
    sample.batches[[a]] <- temp
  }

  sample.batches <- do.list.switch(sample.batches)
  sample.batches

  sample.batches <- do.add.cols(sample.batches, "New", sample.tb, "sample.names")
  sample.batches$RUV_SAMPLES <- paste0(sample.batches$RUV_SAMPLES, "_", sample.batches$Values.RUV_BATCH_NUMBERS)
  sample.batches

  ref.samples <- sample.batches[which(sample.names == ref.samples), "RUV_SAMPLES"][[1]]

  ### Sample and batch preparation

  x <- dat[, c("RUV_SAMPLES", "RUV_CLUSTER_NUMBERS", align.cols), with = FALSE]
  names(x)[c(1:2)] <- c("sample", "cluster")
  x <- as.data.frame(x)

  clust.list <- list(sort(unique(x[["cluster"]])))
  rep.list <- list(ref.samples)

  message(paste0(names(x), " "))
  message(clust.list)
  message(rep.list)

  if (!any(which(unique(dat$RUV_SAMPLES == ref.samples)))) {
    stop("The names of the reference files have been corrupted")
  }

  ### Run RUV

  message("RUV starting")

  res <- run_RUVIII(
    data = x,
    norm_clusters = clust.list,
    k = k,
    rep_samples = rep.list
  )

  res <- as.data.table(res)
  names(res) <- paste0(names(res), append.name)

  ### Finalise and return

  start.dat <- cbind(start.dat, res)

  message("RUV complete")
  return(start.dat)
}

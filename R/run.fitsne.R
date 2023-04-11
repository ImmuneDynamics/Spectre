#' Run FIt-SNE, Fourier Transform TSNE.
#'
#' Implementation of FIt-SNE is available from https://github.com/KlugerLab/FIt-SNE.
#' This function uses \code{\link{fftRtsne}} to run FIt-SNE.
#'
#' @param dat NO DEFAULT. Input data.table or data.frame.
#' @param use.cols NO DEFAULT. Vector of column names or numbers for clustering.
#' @param seed Default = 42. Seed value for reproducibility.
#' @param fitsne.x.name Default = "FItSNE_X". Character. Name of FItSNE x-axis.
#' @param fitsne.y.name Default = "FItSNE_Y". Character. Name of FItSNE y-axis.
#' @param dims Default = 2. Dimensionality of the embedding (reduced data).
#' @param perplexity Default = 30. Perplexity is used to determine the bandwidth of the Gaussian kernel in the input space
#' @param theta Default = 0.5.
#'   For exact t-SNE, set to 0.
#'   If non-zero, then will use either Barnes Hut or FIt-SNE based on nbody_algo.
#'   If Barnes Hut, then this determines the accuracy of BH approximation.
#' @param max_iter Default = 750. Number of iterations of t-SNE to run.
#' @param fft_not_bh Default = TRUE. If theta is nonzero, this determines whether to use FIt-SNE or Barnes Hut approximation.
#' @param ann_not_vptree Default = TRUE. Use vp-trees (as in bhtsne) or approximate nearest neighbors (default).
#'   Set to be TRUE for approximate nearest neighbors.
#' @param exaggeration_factor Default = 12. Coefficient for early exaggeration (>1).
#' @param no_momentum_during_exag Default = FALSE. Set to 0 to use momentum and other optimization tricks.
#'   Can be set to 1 to do plain, vanilla gradient descent (useful for testing large exaggeration coefficients).
#' @param stop_early_exag_iter Default = 250. When to switch off early exaggeration.
#' @param start_late_exag_iter Default = -1. When to start late exaggeration.
#'   Set to -1 by default to not use late exaggeration.
#' @param late_exag_coeff Default = 1. Late exaggeration coefficient.
#'   Set to 1 by default to not use late exaggeration.
#' @param mom_switch_iter Default = 250. Iteration number to switch from momentum to final_momentum.
#' @param momentum Default = 0.5.Initial value of momentum.
#' @param final_momentum Default = 0.8. Value of momentum to use later in the optimisation.
#' @param learning_rate Default = 'auto'. Set to desired learning rate or 'auto', which sets
#'   learning rate to N/exaggeration_factor where N is the sample size, or to
#'   200 if N/exaggeration_factor < 200.
#' @param n_trees Default = 50. When using Annoy, the number of search trees to use.
#' @param search_k Default = -1. When using Annoy, the number of nodes to inspect during
#'   search. Default is -1 which translate to 3*perplexity*n_trees (or K*n_trees when using fixed sigma).
#' @param nterms Default = 3. If using FIt-SNE, this is the number of interpolation points
#'   per sub-interval.
#' @param intervals_per_integer Default = 1. See min_num_intervals.
#' @param min_num_intervals Default = 50. Let maxloc = ceil(max(max(X))) and minloc =
#'   floor(min(min(X))). i.e. the points are in a [minloc]^no_dims by
#'   [maxloc]^no_dims interval/square. The number of intervals in each dimension
#'   is either min_num_intervals or ceil((maxloc -
#'   minloc)/intervals_per_integer), whichever is larger. min_num_intervals must
#'   be an integer >0, and intervals_per_integer must be >0. Defaults are
#'   min_num_intervals=50 and intervals_per_integer = 1.
#' @param sigma Default = -30. Fixed sigma value to use when perplexity==-1.
#' @param K Default = -1. Number of nearest neighbours to get when using fixed sigma.
#' @param initialization Default = 'pca'. pca', 'random', or N x no_dims array to intialize the solution.
#' @param max_step_norm Default = 5. Maximum distance that a point is allowed to move on one
#'   iteration. Larger steps are clipped to this value. This prevents possible
#'   instabilities during gradient descent. Set to -1 to switch it off.
#' @param load_affinities Default = NULL. If 1, input similarities are loaded from a file and
#'   not computed. If 2, input similarities are saved into a file. If 0,
#'   affinities are neither saved nor loaded.
#' @param fast_tsne_path Default = NULL. Path to FItSNE executable.
#' @param nthreads Default = 0. Number of threads to use, set to use all available threads by default.
#' @param perplexity_list Default = NULL. If perplexity==0 then perplexity combination will be
#'   used with values taken from perplexity_list.
#' @param get_costs Default = FALSE. Logical indicating whether the KL-divergence costs computed
#'   every 50 iterations should be returned.
#' @param df Default = 1.0. Positive numeric that controls the degree of freedom of
#'   t-distribution. The actual degree of freedom is 2*df-1. The standard t-SNE
#'   choice of 1 degree of freedom corresponds to df=1. Large df approximates
#'   Gaussian kernel. df<1 corresponds to heavier tails, which can often resolve
#'   substructure in the embedding. See Kobak et al. (2019) for details.
#'
#' @usage
#' run.fitsne(dat, use.cols, seed = 42, fitsne.x.name = "FItSNE_X", fitsne.y.name = "FItSNE_Y",
#'   dims = 2, perplexity = 30, theta = 0.5, max_iter = 750, fft_not_bh = TRUE, ann_not_vptree = TRUE,
#'   stop_early_exag_iter = 250, exaggeration_factor = 12.0, no_momentum_during_exag = FALSE,
#'   start_late_exag_iter = -1, late_exag_coeff = 1.0, mom_switch_iter = 250, momentum = 0.5,
#'   final_momentum = 0.8, learning_rate = 'auto', n_trees = 50, search_k = -1, nterms = 3,
#'   intervals_per_integer = 1, min_num_intervals = 50, K = -1, sigma = -30, initialization = 'pca',
#'   max_step_norm = 5, load_affinities = NULL, fast_tsne_path = NULL, nthreads = 0,
#'   perplexity_list = NULL, get_costs = FALSE,  df = 1.0)
#'
#' @examples
#' dat <- demo.clustered
#' dat.sub <- do.subsample(dat, 30000)
#' use.cols <- names(dat)[12:19]
#' dat.reduced <- run.fitsne(dat = dat.sub, use.cols = use.cols)
#'
#' @author
#' Givanna Putri
#' 
#' @import data.table
#'
#' @export
run.fitsne <- function(dat,
                       use.cols,
                       seed = 42,
                       fitsne.x.name = "FItSNE_X",
                       fitsne.y.name = "FItSNE_Y",
                       dims = 2,
                       perplexity = 30,
                       theta = 0.5,
                       max_iter = 750,
                       fft_not_bh = TRUE,
                       ann_not_vptree = TRUE,
                       stop_early_exag_iter = 250,
                       exaggeration_factor = 12.0,
                       no_momentum_during_exag = FALSE,
                       start_late_exag_iter = -1,
                       late_exag_coeff = 1.0,
                       mom_switch_iter = 250,
                       momentum = 0.5,
                       final_momentum = 0.8,
                       learning_rate = "auto",
                       n_trees = 50,
                       search_k = -1,
                       nterms = 3,
                       intervals_per_integer = 1,
                       min_num_intervals = 50,
                       K = -1,
                       sigma = -30,
                       initialization = "pca",
                       max_step_norm = 5,
                       load_affinities = NULL,
                       fast_tsne_path = NULL,
                       nthreads = 0,
                       perplexity_list = NULL,
                       get_costs = FALSE,
                       df = 1.0) {

  if (is.null(fast_tsne_path)) {
    # use pre-compiled one given in Spectre
    if (.Platform$OS.type == "unix") {
      fast_tsne_path <- file.path(.libPaths()[1], "Spectre", "FItSNE_executable", "FItSNE_Unix_1_2_1")
    } else if (.Platform$OS.type == "windows") {
      fast_tsne_path <- file.path(.libPaths()[1], "Spectre", "FItSNE_executable", "FItSNE-Windows-1.2.1", "FItSNE.exe")
    }
  }

  fitsne_out <- fftRtsne(
    X = as.matrix(dat[, use.cols, with = FALSE]),
    rand_seed = seed,
    dims = dims,
    perplexity = perplexity,
    theta = theta,
    max_iter = max_iter,
    fft_not_bh = fft_not_bh,
    ann_not_vptree = ann_not_vptree,
    stop_early_exag_iter = stop_early_exag_iter,
    exaggeration_factor = exaggeration_factor,
    no_momentum_during_exag = no_momentum_during_exag,
    start_late_exag_iter = start_late_exag_iter,
    late_exag_coeff = late_exag_coeff,
    mom_switch_iter = mom_switch_iter,
    momentum = momentum,
    final_momentum = final_momentum,
    learning_rate = learning_rate,
    n_trees = n_trees,
    search_k = search_k,
    nterms = nterms,
    intervals_per_integer = intervals_per_integer,
    min_num_intervals = min_num_intervals,
    K = K,
    sigma = sigma,
    initialization = initialization,
    max_step_norm = max_step_norm,
    load_affinities = load_affinities,
    fast_tsne_path = fast_tsne_path,
    nthreads = nthreads,
    perplexity_list = perplexity_list,
    get_costs = get_costs,
    df = df
  )
  dat_bk <- copy(dat)
  dat_bk[[fitsne.x.name]] <- fitsne_out[, 1]
  dat_bk[[fitsne.y.name]] <- fitsne_out[, 2]

  return(dat_bk)
}

#' FIt-SNE Based on Kluger Lab FIt-SNE
#'
#' Modified version of \url{https://github.com/KlugerLab/FIt-SNE/}
#' Please do not directly call this function as it requires compiled FIt-SNE code to run.
#' If you want to run FIt-SNE, please have a look at \code{\link{run.fitsne}} function.
#'
#' @seealso \code{\link{run.fitsne}}
#'
#' @references Linderman, G., Rachh, M., Hoskins, J., Steinerberger, S.,
#'   Kluger., Y. (2019). Fast interpolation-based t-SNE for improved
#'   visualization of single-cell RNA-seq data. Nature Methods.
#'   \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6402590/}.
#' 
#' @import rsvd irlba
#' @noRd
#' 
fftRtsne <- function(X,
                     dims = 2,
                     perplexity = 30,
                     theta = 0.5,
                     max_iter = 750,
                     fft_not_bh = TRUE,
                     ann_not_vptree = TRUE,
                     stop_early_exag_iter = 250,
                     exaggeration_factor = 12.0,
                     no_momentum_during_exag = FALSE,
                     start_late_exag_iter = -1,
                     late_exag_coeff = 1.0,
                     mom_switch_iter = 250,
                     momentum = 0.5,
                     final_momentum = 0.8,
                     learning_rate = "auto",
                     n_trees = 50,
                     search_k = -1,
                     rand_seed = -1,
                     nterms = 3,
                     intervals_per_integer = 1,
                     min_num_intervals = 50,
                     K = -1,
                     sigma = -30,
                     initialization = "pca",
                     max_step_norm = 5,
                     data_path = NULL,
                     result_path = NULL,
                     load_affinities = NULL,
                     fast_tsne_path = NULL,
                     nthreads = 0,
                     perplexity_list = NULL,
                     get_costs = FALSE,
                     df = 1.0) {
    version_number <- "1.2.1"
    
    if (is.null(data_path)) {
        data_path <- tempfile(pattern = 'fftRtsne_data_', fileext = '.dat')
    }
    if (is.null(result_path)) {
        result_path <- tempfile(pattern = 'fftRtsne_result_', fileext = '.dat')
    }
    
    if (is.null(fast_tsne_path)) {
        fast_tsne_path <- system2('which', 'fast_tsne', stdout = TRUE)
    }
    fast_tsne_path <- normalizePath(fast_tsne_path)
    if (!file_test("-x", fast_tsne_path)) {
        stop(fast_tsne_path, " does not exist or is not executable; check your fast_tsne_path parameter")
    }
    
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    
    if (!is.numeric(theta) || (theta < 0.0) || (theta > 1.0) ) { stop("Incorrect theta.")}
    if (nrow(X) - 1 < 3 * perplexity) { stop("Perplexity is too large.")}
    if (!is.matrix(X)) { stop("Input X is not a matrix")}
    if (!(max_iter > 0)) { stop("Incorrect number of iterations.")}
    if (!is.wholenumber(stop_early_exag_iter) || stop_early_exag_iter < 0) { stop("stop_early_exag_iter should be a positive integer")}
    if (!is.numeric(exaggeration_factor)) { stop("exaggeration_factor should be numeric")}
    if (!is.numeric(df)) { stop("df should be numeric")}
    if (!is.wholenumber(dims) || dims <= 0) { stop("Incorrect dimensionality.")}
    if (search_k == -1) {
        if (perplexity > 0) {
            search_k <- n_trees * perplexity * 3
        } else if (perplexity == 0) {
            search_k <- n_trees * max(perplexity_list) * 3
        } else { 
            search_k <- n_trees * K
        }
    }
    
    if (is.character(learning_rate) && learning_rate =='auto') {
        learning_rate = max(200, nrow(X)/exaggeration_factor)
    }
    if (is.character(start_late_exag_iter) && start_late_exag_iter =='auto') {
        if (late_exag_coeff > 0) {
            start_late_exag_iter = stop_early_exag_iter
        }else {
            start_late_exag_iter = -1
        }
    }
    
    if (is.character(initialization) && initialization =='pca') {
        if (rand_seed != -1)  {
            set.seed(rand_seed)
        }
        if ("rsvd" %in% utils::installed.packages()) {
            message('Using rsvd() to compute the top PCs for initialization.')
            X_c <- scale(X, center=T, scale=F)
            rsvd_out <- rsvd::rsvd(X_c, k=dims)
            X_top_pcs <- rsvd_out$u %*% diag(rsvd_out$d, nrow=dims)
        }else if("irlba" %in% utils::installed.packages()) { 
            message('Using irlba() to compute the top PCs for initialization.')
            X_colmeans <- colMeans(X)
            irlba_out <- irlba::irlba(X,nv=dims, center=X_colmeans)
            X_top_pcs <- irlba_out$u %*% diag(irlba_out$d, nrow=dims)
        }else{
            stop("By default, FIt-SNE initializes the embedding with the
                     top PCs. We use either rsvd or irlba for fast computation.
                     To use this functionality, please install the rsvd package
                     with install.packages('rsvd') or the irlba package with
                     install.packages('ilrba').  Otherwise, set initialization
                     to NULL for random initialization, or any N by dims matrix
                     for custom initialization.")
        }
        initialization <- 0.0001*(X_top_pcs/sd(X_top_pcs[,1])) 
        
    }else if (is.character(initialization) && initialization == 'random'){
        message('Random initialization')
        initialization = NULL
    }
    
    if (fft_not_bh) {
        nbody_algo <- 2
    } else {
        nbody_algo <- 1
    }
    
    if (is.null(load_affinities)) {
        load_affinities <- 0
    } else {
        if (load_affinities == 'load') {
            load_affinities <- 1
        } else if (load_affinities == 'save') {
            load_affinities <- 2
        } else {
            load_affinities <- 0
        }
    }
    
    if (ann_not_vptree) {
        knn_algo <- 1
    } else {
        knn_algo <- 2
    }
    tX <- as.numeric(t(X))
    
    f <- file(data_path, "wb")
    n <- nrow(X)
    D <- ncol(X)
    writeBin(as.integer(n), f, size = 4)
    writeBin(as.integer(D), f, size = 4)
    writeBin(as.numeric(theta), f, size = 8) #theta
    writeBin(as.numeric(perplexity), f, size = 8)
    
    if (perplexity == 0) {
        writeBin(as.integer(length(perplexity_list)), f, size = 4)
        writeBin(perplexity_list, f) 
    }
    
    writeBin(as.integer(dims), f, size = 4)
    writeBin(as.integer(max_iter), f, size = 4)
    writeBin(as.integer(stop_early_exag_iter), f, size = 4)
    writeBin(as.integer(mom_switch_iter), f, size = 4)
    writeBin(as.numeric(momentum), f, size = 8)
    writeBin(as.numeric(final_momentum), f, size = 8)
    writeBin(as.numeric(learning_rate), f, size = 8)
    writeBin(as.numeric(max_step_norm), f, size = 8)
    writeBin(as.integer(K), f, size = 4) #K
    writeBin(as.numeric(sigma), f, size = 8) #sigma
    writeBin(as.integer(nbody_algo), f, size = 4)  #not barnes hut
    writeBin(as.integer(knn_algo), f, size = 4) 
    writeBin(as.numeric(exaggeration_factor), f, size = 8) #compexag
    writeBin(as.integer(no_momentum_during_exag), f, size = 4) 
    writeBin(as.integer(n_trees), f, size = 4) 
    writeBin(as.integer(search_k), f, size = 4) 
    writeBin(as.integer(start_late_exag_iter), f, size = 4) 
    writeBin(as.numeric(late_exag_coeff), f, size = 8) 
    
    writeBin(as.integer(nterms), f, size = 4) 
    writeBin(as.numeric(intervals_per_integer), f, size = 8) 
    writeBin(as.integer(min_num_intervals), f, size = 4) 
    writeBin(tX, f) 
    writeBin(as.integer(rand_seed), f, size = 4) 
    writeBin(as.numeric(df), f, size = 8)
    writeBin(as.integer(load_affinities), f, size = 4) 
    if (!is.null(initialization) ) { writeBin( c(t(initialization)), f) }		
    close(f) 
    
    flag <- system2(command = fast_tsne_path, 
                    args = c(version_number, data_path, result_path, nthreads))
    if (flag != 0) {
        stop('tsne call failed')
    }
    f <- file(result_path, "rb")
    n <- readBin(f, integer(), n = 1, size = 4)
    d <- readBin(f, integer(), n = 1, size = 4)
    Y <- readBin(f, numeric(), n = n * d)
    Y <- t(matrix(Y, nrow = d))
    if (get_costs) {
        readBin(f, integer(), n = 1, size = 4)
        costs <- readBin(f, numeric(), n = max_iter, size = 8)
        Yout <- list(Y = Y, costs = costs)
    } else {
        Yout <- Y
    }
    close(f)
    file.remove(data_path)
    file.remove(result_path)
    Yout
}

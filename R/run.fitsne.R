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
  if (!is.element("rsvd", installed.packages()[, 1])) stop("rsvd is required but not installed")
  if (!is.element("irlba", installed.packages()[, 1])) stop("irlba is required but not installed")
  if (!is.element("data.table", installed.packages()[, 1])) stop("data.table is required but not installed")

  require("rsvd")
  require("irlba")

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

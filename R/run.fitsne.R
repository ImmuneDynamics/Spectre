#' Run FIt-SNE, Fourier Transform TSNE.
#' 
#' Implementation of FIt-SNE is available from https://github.com/KlugerLab/FIt-SNE.
#' As it is currently not possible to install FIt-SNE directly from github using devtools,
#' we made a verbatim copy of KlugerLab's fftRtsne, and source it in this method.
#' 
#' 
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
                       learning_rate = 'auto',
                       n_trees = 50, 
                       search_k = -1, 
                       nterms = 3, 
                       intervals_per_integer = 1, 
                       min_num_intervals = 50, 
                       K = -1, 
                       sigma = -30, 
                       initialization = 'pca',
                       max_step_norm = 5,
                       load_affinities = NULL, 
                       fast_tsne_path = NULL,
                       nthreads = 0, 
                       perplexity_list = NULL, 
                       get_costs = FALSE, 
                       df = 1.0) {
    
    if(!is.element('rsvd', installed.packages()[,1])) stop('rsvd is required but not installed')
    if(!is.element('ilrba', installed.packages()[,1])) stop('ilrba is required but not installed')
    if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')
    
    require('rsvd')
    require('ilrba')
    
    if (is.null(fast_tsne_path)) {
        # use pre-compiled one given in Spectre
        if (.Platform$OS.type == "unix") {
            fast_tsne_path <- file.path(.libPaths(), "Spectre", "bin", "fast_tsne")
        }
        # need to change the following to libpath
        # and get the exe version
        # else {
        #     fast_tsne_path <- file.path(FAST_TSNE_SCRIPT_DIR, "bin", "FItSNE.exe")
        # }
    }
    
    fitsne_out <- fftRtsne(X = as.matrix(dat[, use.cols, with = FALSE]),
                           rand_seed = seed,
                           fast_tsne_path = fast_tsne_path)
    dat_bk <- copy(dat)
    dat_bk[[fitsne.x.name]] <- fitsne_out[,1]
    dat_bk[[fitsne.y.name]] <- fitsne_out[,2]
    
    return(dat_bk)
    
}
#' Title
#'
#' @param dat 
#' @param markers 
#' @param sample_source_col 
#' @param batch_id_col 
#'
#' @return
#' @export
#'
#' @importFrom graphics hist
#' @importFrom emdist emd2d
#' @import data.table
#' 
do.calculate.earth.mover.dist <- function(dat, markers, sample_source_col, batch_id_col) {
    
    sample_sources <- unique(dat[[sample_source_col]])
    batches <- unique(dat[[batch_id_col]])
    
    emd_dist <- lapply(sample_sources, function(sample_source) {
        
        emd_per_marker <- lapply(markers, function(marker) {
            
            exp <- dat[get(sample_source_col) == sample_source, c(batch_id_col, marker), with = FALSE]
            
            hist_bins <- lapply(batches, function(bat_id) {
                exp_bat <- exp[get(batch_id_col) == bat_id, ][[marker]]
                exp_bat <- as.matrix(hist(exp_bat, breaks = seq(-100, 100, by = 0.1), plot = FALSE)$counts)
                return(exp_bat)
            })
            emd2d(hist_bins[[1]], hist_bins[[2]])
            
        })
        emd_per_marker <- data.table(
            marker = markers,
            emd = unlist(emd_per_marker),
            sample_source = sample_source
        )
        setnames(emd_per_marker, "sample_source", "sample_source_col")
        return(emd_per_marker)
        
    })
    emd_dist <- rbindlist(emd_dist)
    return(emd_dist)
}

#' Draw MDS plot
#' 
#' Experimental feature. Still very much work in progress.
#'
#' @param dat 
#' @param sample_col 
#' @param markers 
#' @param colour_by 
#'
#' @return
#' @export
#'
#' @importFrom limma plotMDS
#' @import data.table
#' @import ggplot2 
make.mds.plot <- function(dat,
                          sample_col,
                          markers,
                          colour_by) {
    pseudobulk <- dat[, lapply(.SD, mean), by = sample_col, .SDcols = markers]
    pseudobulk_transposed <- t(pseudobulk[, markers, with = FALSE])
    colnames(pseudobulk_transposed) <- pseudobulk[[sample_col]]
    
    mds <- plotMDS(pseudobulk_transposed)
    
    # To get unique combination of sample and colour by
    sample_colour_by_combo <- unique(dat[, c(sample_col, colour_by), with = FALSE])
    
    mds_df <- merge.data.table(
        x = data.frame(
            dim1 = mds$x,
            dim2 = mds$y,
            sample_id = pseudobulk[[sample_col]]
        ),
        y = sample_colour_by_combo,
        by = sample_col
    )
    
    # convert colour by column to factor
    mds_df[[colour_by]] <- factor(mds_df[[colour_by]])
    
    
    plt <- ggplot(mds_df, aes(x = dim1, y = dim2, color = !! sym(colour_by))) +
        geom_point() +
        geom_label_repel(aes(label = !! sym(sample_col)), show.legend = FALSE, max.overlaps = 20) +
        theme_bw() +
        labs(
            x = paste0("MDS1 (", round( 100 * mds$var.explained[1]), "%)"), 
            y = paste0("MDS2 (", round( 100 * mds$var.explained[2]), "%)"),
            color = colour_by,
            title = "MDS plot of pseudobulk protein expression of samples")
    
    return(plt)
    
}
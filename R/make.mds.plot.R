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
                          colour_by,
                          font_size = 4) {
    agg_dat <- dat[, lapply(.SD, mean), by = sample_col, .SDcols = markers]
    agg_dat_transposed <- t(agg_dat[, markers, with = FALSE])
    colnames(agg_dat_transposed) <- agg_dat[[sample_col]]
    
    mds <- limma::plotMDS(agg_dat_transposed, plot = FALSE)
    
    # To get unique combination of sample and colour by
    sample_colour_by_combo <- unique(dat[, c(sample_col, colour_by), with = FALSE])
    
    # Create data frame
    mds.df <- data.frame(
      dim1 = mds$x,
      dim2 = mds$y,
      placeholder.name = agg_dat[[sample_col]]
    )
    
    # Replace column name (important for merging with original data)
    colnames(mds.df) <- sub("placeholder.name", sample_col, colnames(mds.df))
    
    mds_df <- merge.data.table(
        x = mds.df,
        y = sample_colour_by_combo,
        by = sample_col
    )
    
    # convert colour by column to factor
    mds_df[[colour_by]] <- factor(mds_df[[colour_by]])
    
    
    plt <- ggplot(mds_df, aes(x = dim1, y = dim2, color = !! sym(colour_by))) +
        geom_point() +
        ggrepel::geom_label_repel(aes(label = !! sym(sample_col)), 
                         show.legend = FALSE, max.overlaps = 50, size = font_size, 
                         box.padding = unit(0.1, "lines")) +
    theme_bw() +
    labs(
        x = paste0("MDS1 (", round( 100 * mds$var.explained[1]), "%)"), 
        y = paste0("MDS2 (", round( 100 * mds$var.explained[2]), "%)"),
        color = colour_by,
        title = "MDS plot of mean protein expression of samples")

return(plt)

}
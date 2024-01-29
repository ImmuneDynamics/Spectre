#' Make Bubble Plot
#' 
#' Dot (or bubble) plot which visualises the marker expression changes and the 
#' percentage of cells captured across different clusters.
#' Experimental feature.
#' 
#' @param dt Data.table object containing the data to plot
#' @param markers_to_plot Vector of markers to plot
#' @param cluster Name of the column in dt which specify the clusters
#' @col.scheme Colour scheme to use. Must be something from 
#' the RColorBrewer or the viridis package.
#' Default to 'RdYlBu' in the RColorBrewer package
#' 
#' @return ggplot object containing the plot
#' 
#' @import data.table
#' @import ggplot2
#' @import viridis
#' @importFrom RColorBrewer brewer.pal.info
#' 
#' @author Givanna Putri
#' @export

make.bubble.plot <- function(
        dt, 
        markers_to_plot, 
        cluster, 
        col.scheme='RdYlBu'
) {
    # require
    # ggplot2, rcolorbrewer, viridis
    
    # Calculate clusters' markers' mean expression and proportion of cells captured.
    avg_exp <- dt[, lapply(.SD, mean), .SDcols = markers_to_plot, by = cluster]
    cluster_prop <- dt[, .(prop = .N * 100/nrow(dt)), by = cluster]
    
    dt_for_plt <- melt(avg_exp, id.vars=c(cluster), variable.name='marker')
    dt_for_plt <- merge.data.table(dt_for_plt, cluster_prop, by=cluster)
    dt_for_plt[[cluster]] <- factor(dt_for_plt[[cluster]])
    
    
    plt <- ggplot(dt_for_plt, aes(x = marker, y = .data[[cluster]])) +
        geom_point(aes(color = value, size = prop)) +
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        labs(x = 'Marker', size = 'Percent of Cells', 
             color = 'Mean Expression')
    
    # Work out what colour scheme was selected.
    # With RColorBrewer, it is easy as we can just check whether
    # the colour exists using brewer.pal.info
    
    if (!is.na(brewer.pal.info[col.scheme, "maxcolors"])) {
        plt <- plt + scale_color_distiller(palette = col.scheme)
    } else {
        warning(paste("Finding col.scheme", col.scheme, "in viridis package."))
        plt <- plt + scale_color_viridis(option = col.scheme)
    }
    return(plt)
}

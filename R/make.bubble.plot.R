#' Make Bubble Plot
#' 
#' Dot (or bubble) plot which visualises the marker expression changes and the 
#' percentage of cells captured across different clusters.
#' Experimental feature.
#' Use at your own risk.
#' 
#' @param dt Data.table object containing the data to plot
#' @param markers_to_plot Vector of markers to plot
#' @param cluster Name of the column in dt which specify the clusters
#' @col.scheme Colour scheme to use. Must be something from the viridis package.
#' Default to 'viridis'
#' 
#' @return ggplot object containing the plot
#' 
#' @export make.bubble.plot

make.bubble.plot <- function(dt, markers_to_plot, cluster, col.scheme='viridis') {
    
    require(data.table)
    require(ggplot2)
    require(viridis)
    
    avg_exp <- dt[, lapply(.SD, mean), .SDcols = markers_to_plot, by = cluster]
    cluster_prop <- dt[, .(prop = .N * 100/nrow(dt)), by = cluster]
    
    dt_for_plt <- melt(avg_exp, id.vars=c(cluster), variable.name='marker')
    dt_for_plt <- merge.data.table(dt_for_plt, cluster_prop, by=cluster)
    dt_for_plt[[cluster]] <- factor(dt_for_plt[[cluster]])
    
    plt <- ggplot(dt_for_plt, aes_string(x = "marker", y = cluster)) +
        geom_point(aes(color = value, size = prop)) +
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        labs(x = 'Marker', size = 'Percent of Cells', 
             color = 'Average Expression') +
        scale_color_viridis(option=col.scheme)
    return(plt)
}
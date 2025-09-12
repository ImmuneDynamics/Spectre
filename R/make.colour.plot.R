#' Make colour plot
#' 
#' Create a dot plot (e.g., UMAP or tSNE) of cells, 
#' coloured by a continuous variable (e.g., marker expression) or a factor 
#' (e.g., cluster, group).
#' 
#' @param dat A data.table containing the data to plot.
#' @param x.axis Character. Column name for the x-axis.
#' @param y.axis Character. Column name for the y-axis.
#' @param col.axis Character or NULL. Column name for colouring points. 
#' If NULL, points are coloured by density.
#' @param col.type Character. "continuous" (default) or "factor". 
#' Determines how \code{col.axis} is interpreted.
#' @param add.label Logical. If TRUE and \code{col.type = "factor"}, 
#' adds labels at the centroid of each group.
#' @param hex Logical. If TRUE, uses hex binning (only for continuous colour plots).
#' @param hex.bins Integer. Number of hex bins if \code{hex = TRUE}.
#' @param colours Character. Colour scheme for continuous plots. 
#' Options available are: "jet", all options in \code{RColorBrewer::brewer.pal.info},
#' and all options in viridis pallete. 
#' Default is "spectral".
#' @param col.min.threshold Numeric. Minimum quantile for colour scale (continuous).
#' @param col.max.threshold Numeric. Maximum quantile for colour scale (continuous).
#' @param align.xy.by data.table. Data to use for aligning x/y axis limits.
#' @param align.col.by data.table. Data to use for aligning colour scale limits.
#' @param regression.line Character or NULL. If not NULL, 
#' adds a regression line ("lm", "loess", etc.).
#' @param title Character or NULL. Plot title. Defaults to \code{col.axis}.
#' @param filename Character or NULL. File name for saving the plot. 
#' If NULL, a name is generated automatically.
#' @param dot.size Numeric. Size of points.
#' @param plot.width Numeric. Width of saved plot (in inches).
#' @param plot.height Numeric. Height of saved plot (in inches).
#' @param nudge_x, nudge_y Numeric. Amount to nudge centroid labels 
#' (if \code{add.label = TRUE}).
#' @param square Logical. If TRUE, enforces a square aspect ratio.
#' @param legend.loc Character. 
#' Legend position: "right" (default), "bottom", "top", "left", or "none".
#' @param save.to.disk Logical. If TRUE (default), saves the plot to disk. If FALSE, only displays the plot.
#' @param path Character. Directory to save the plot.
#' @param blank.axis Logical. If TRUE, produces a minimalist plot with no axis lines or labels.
#'
#' @return A ggplot2 object representing the plot.

#' @usage make.colour.plot(dat, x.axis, y.axis)
#'
#' @examples
#' Spectre::make.colour.plot(
#'     dat = Spectre::demo.clustered,
#'     x.axis = "UMAP_X",
#'     y.axis = "UMAP_Y",
#'     col.axis = "CD4_asinh"
#' )
#'
#' @author
#' Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Givanna Putri
#'
#' @import data.table
#' @import ggplot2
#' @import scales
#' @import colorRamps
#' @import ggthemes
#' @import RColorBrewer
#' @import ggrepel
#'
#' @export
#' 
make.colour.plot <- function(
    dat,
    x.axis,
    y.axis,
    col.axis = NULL,
    col.type = c("continuous", "factor"),
    add.label = FALSE,
    hex = FALSE,
    hex.bins = 30,
    colours = "spectral",
    col.min.threshold = 0.01,
    col.max.threshold = 0.995,
    align.xy.by = dat,
    align.col.by = dat,
    regression.line = NULL,
    title = col.axis,
    filename = NULL,
    dot.size = 1,
    plot.width = 9,
    plot.height = 7,
    nudge_x = 0.5,
    nudge_y = 0.5,
    square = TRUE,
    legend.loc = c("right", "bottom", "top", "left", "none"),
    save.to.disk = TRUE,
    path = getwd(),
    blank.axis = FALSE
) {
    
    # For testing
    # dat <- Spectre::demo.clustered
    # x.axis <- 'UMAP_X'
    # y.axis <- 'UMAP_Y'
    # col.axis <- 'Population'
    # col.type = "continuous" # can be "continuous" or "factor"
    # add.label = FALSE # only works for 'factor'
    # hex = FALSE
    # hex.bins = 30
    # colours = "spectral" # can be spectral, jet, etc      # only works for continuous
    # col.min.threshold = 0.01
    # col.max.threshold = 0.995
    # align.xy.by = dat
    # align.col.by = dat
    # regression.line = NULL # "lm" # "loess"
    # title = col.axis
    # filename = NULL
    # dot.size = 1
    # plot.width = 9
    # plot.height = 7
    # nudge_x = 0.5
    # nudge_y = 0.5
    # square = TRUE
    # legend.loc = NULL # 'right' and 'bottom'
    # save.to.disk = TRUE
    # path = getwd()
    # blank.axis = FALSE

    col.type <- tryCatch(
        match.arg(col.type),
        error = function(e) {
            stop("Invalid value for 'col.type': must be either 'continuous' or 'factor'.", call. = FALSE)
        }
    )

    # some checks

    if (hex == TRUE) {
        if (is.null(col.axis)) {
            message("Note: hex bins do not currently work for density plots, only for colour plots when col.axis is specified and can be plotted as a continuous numeric variable")
        } else {
            if (!is.numeric(dat[[col.axis]])) { # tests to see if there are any non numeric values
                stop("Sorry, hex bins only work when col.axis is specified, and can be plotted as a continuous numeric variable")
            }
        }
    }

    if (!is.null(col.axis)) {
        if (col.type == "continuous") {
            if (!is.numeric(dat[[col.axis]])) {
                message("Non-numeric values detected in col.axis -- changing col.type to 'factor'")
                col.type <- "factor"
            }
        }

        if (col.type == "factor") {
            if (length(unique(as.factor(dat[[col.axis]]))) > 200) {
                message("Over 200 different factors detected, using continuous scale instead of a factor scale")
                col.type <- "continuous"
            }
        }
    }

    ### Initialise plot

    if (!is.null(col.axis)) {
        if (col.type == "continuous") {
            p <- .make_continuous_scatter_plot(dat, x.axis, y.axis, col.axis, colours, dot.size, hex, hex.bins, align.col.by, col.min.threshold, col.max.threshold)
        } else if (col.type == "factor") {
            p <- .make_factor_scatter_plot(dat, x.axis, y.axis, col.axis, align.col.by, dot.size)
        }
    } else {
        p <- .make_density_plot(dat, x.axis, y.axis, dot.size, colours)
    }

    ### Regression line
    if (!is.null(regression.line)) {
        p <- p + geom_smooth(method = regression.line)
    }

    ### Add title

    if (is.null(title)) {
        title <- "Density"
    }

    p <- p + ggtitle(title)

    ### Set up axis
    ### Define limits for x and y axis
    axis.lims <- .get_axis_limits(dat, x.axis, y.axis, align.xy.by)
    p <- p + scale_x_continuous(breaks = scales::pretty_breaks(n = 8), name = x.axis, limits = c(axis.lims$Xmin, axis.lims$Xmax))
    p <- p + scale_y_continuous(breaks = scales::pretty_breaks(n = 8), name = y.axis, limits = c(axis.lims$Ymin, axis.lims$Ymax))

    ### Set up themes etc
    if (blank.axis) {
        p <- p + theme(
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank(),
        )
    } else {
        p <- p + theme(
            panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
            axis.title.x = element_text(color = "Black", size = 28),
            axis.title.y = element_text(color = "Black", size = 28),
            axis.text.x = element_text(color = "Black", size = 24),
            axis.text.y = element_text(color = "Black", size = 24),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
            plot.title = element_text(color = "Black", face = "bold", size = 32, hjust = 0)
        )
    }
    

    if (square == TRUE) {
        p <- p + theme(aspect.ratio = 1)
    }

    ### Setup legend
    legend.loc <- tryCatch(
        match.arg(legend.loc),
        error = function(e) {
            stop("Invalid value for 'legend.loc': must be either 'top', 'bottom', 'left', 'right', or 'none'.", call. = FALSE)
        }
    )
    if (legend.loc %in% c("top", "bottom")) {
        legend.direction <- "horizontal"
    } else {
        legend.direction <- "vertical"
    }
    p <- p + theme(
        legend.direction = legend.direction,
        legend.position = legend.loc,
        legend.text = element_text(size = 18),
        legend.title = element_blank()
    )

    # add centroid
    if(add.label && col.type == "factor") {
        centroids <- .calc_centroids_dt(dat, x.axis, y.axis, col.axis)

        p <- p + ggplot2::geom_point(
            data = centroids, 
            aes(x = centroidX, y = centroidY), 
            color = "black", size = 2
        )
        p <- p + ggrepel::geom_label_repel(
            data = centroids,
            aes(x = centroidX, y = centroidY, label = centroidCol),
            color = "black",
            fontface = "bold",
            nudge_x = nudge_x,
            nudge_y = nudge_y,
            alpha = 0.8
        )
        p <- p + ggplot2::guides(alpha = "none")
        
    }

    ### Save plot
    if (save.to.disk) {
        if (is.null(filename)) {
            if (!is.null(col.axis)) {
                if (col.type == "continuous") {
                    lb <- "Colour"
                }

                if (col.type == "factor") {
                    lb <- "Factor"
                }
            } else {
                lb <- "Density"
            }
            filename <- paste0(lb, " plot - ", title, " - plotted on ", x.axis, " by ", y.axis, ".png")
        }

        ggplot2::ggsave(
            filename = filename,
            plot = p,
            path = path,
            width = plot.width,
            height = plot.height,
            limitsize = FALSE
        )
    } else {
        print(p)
    }

    return(p)
}


#' Generate a colour scheme
#'
#' Internal function creates a colour scheme based on the provided colours.
#'
#' @param colour_scheme A character string specifying the colour scheme to use.
#'
#' @return A colour scheme object or vector, depending on implementation.
#' @keywords internal
#' @noMd 
#'
.get_colour_scheme <- function(colour_scheme) {

    brewer_opts <- RColorBrewer::brewer.pal.info
    viridis_opts <- c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo")

    # put spectral here just for backward compatibility
    available_colours <- c(viridis_opts, rownames(brewer_opts), "jet", "spectral")
    if (!colour_scheme %in% available_colours) {
        warning(
            sprintf(
                "Invalid colour_scheme '%s'. Valid options are: %s. Defaulting to 'spectral'",
                colour_scheme,
                paste(
                    available_colours,
                    collapse = ", "
                )
            ),
            call. = FALSE
        )
        colour_scheme <- "spectral"
    }

    if (colour_scheme %in% viridis_opts) {
        res <- colorRampPalette(c(viridis_pal(option = colour_scheme)(50)))
    } else if (colour_scheme == "jet") {
        res <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    } else if (tolower(colour_scheme) == "spectral") {
        spectral.list <- colorRampPalette(brewer.pal(11, "Spectral"))(50)
        spectral.list <- rev(spectral.list)
        res <- colorRampPalette(c(spectral.list))
    } else if (colour_scheme %in% rownames(brewer_opts)) {
        n_max_colour <- brewer_opts[colour_scheme, "maxcolors"]
        colour.list <- colorRampPalette(RColorBrewer::brewer.pal(n_max_colour, colour_scheme))(50)
        res <- colorRampPalette(c(colour.list))
    }
    res
}

#' Calculate Axis Limits for Plotting
#'
#' Computes the minimum and maximum values for the x and y axes from the provided data,
#' with optional alignment of axes by a specified variable.
#'
#' @param dat A data frame containing the data to be plotted.
#' @param x.axis A string specifying the column name in \code{dat} to use for the x-axis.
#' @param y.axis A string specifying the column name in \code{dat} to use for the y-axis.
#' @param align.xy.by Optional. A string specifying a column name in \code{dat} by which to align the x and y axes.
#'
#' @return A list containing the calculated limits for the x and y axes.
#'
#' @keywords internal
#' @noMd 
#' @import data.table
#' 
.get_axis_limits <- function(dat, x.axis, y.axis, align.xy.by) {
    # Returns a named list with Xmin, Xmax, Ymin, Ymax
    if (is.null(align.xy.by)) {
        Xmax <- max(dat[[x.axis]], na.rm = TRUE)
        Xmin <- min(dat[[x.axis]], na.rm = TRUE)
        Ymax <- max(dat[[y.axis]], na.rm = TRUE)
        Ymin <- min(dat[[y.axis]], na.rm = TRUE)
    } else {
        Xmax <- max(align.xy.by[[x.axis]], na.rm = TRUE)
        Xmin <- min(align.xy.by[[x.axis]], na.rm = TRUE)
        Ymax <- max(align.xy.by[[y.axis]], na.rm = TRUE)
        Ymin <- min(align.xy.by[[y.axis]], na.rm = TRUE)
    }
    list(Xmin = Xmin, Xmax = Xmax, Ymin = Ymin, Ymax = Ymax)
}

#' Create a density scatter plot with customizable dot size and color scheme
#'
#' This internal function generates a density scatter plot from the provided data frame,
#' plotting the specified variables on the x and y axes. The appearance of the plot can be
#' customized by adjusting the dot size and selecting a color palette.
#'
#' @param dat A data frame containing the data to be plotted.
#' @param x.axis A string specifying the column name in \code{dat} to be used for the x-axis.
#' @param y.axis A string specifying the column name in \code{dat} to be used for the y-axis.
#' @param dot.size Numeric value indicating the size of the dots in the plot. Default is 1.
#' @param colours A string specifying the color palette to use. Default is "spectral".
#'
#' @return A ggplot2 object representing the density scatter plot.
#' @keywords internal
#' @noMd 
#' 
#' @import ggplot2
#' @import ggpointdensity
#' @import viridis
#' @import RColorBrewer
#' 
.make_density_plot <- function(dat, x.axis, y.axis, dot.size, colours) {
    p <- ggplot2::ggplot(
        data = dat,
        aes(
            x = .data[[x.axis]],
            y = .data[[y.axis]]
        )
    ) +
        ggpointdensity::geom_pointdensity(size = dot.size)

    viridis_opts <- c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo")

    if (colours %in% viridis_opts) {
        p <- p + viridis::scale_colour_viridis(option = colours)
    } else if (colours == "jet") {
        p <- p + ggplot2::scale_colour_gradientn(colours = c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    } else if (colours == "spectral") {
        p <- p + ggplot2::scale_colour_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(50)))
    } else if (colours == "BuPu") {
        colour.list <- colorRampPalette(RColorBrewer::brewer.pal(9, "BuPu"))(31)
        p <- p + ggplot2::scale_colour_gradientn(colours = colour.list)
    }
    return(p)
}

#' Internal: Create a factor colour scatter plot
#'
#' Generates a scatter plot for factor (categorical) colouring.
#'
#' @param dat Data frame or data.table.
#' @param x.axis Character, column name for x axis.
#' @param y.axis Character, column name for y axis.
#' @param col.axis Character, column name for colour (factor).
#' @param dot.size Numeric, size of dots.
#'
#' @return ggplot2 object
#' @keywords internal
#' @noMd 
#' 
#' @import ggplot2
#' @import data.table
#' 
.make_factor_scatter_plot <- function(
    dat, x.axis, y.axis, col.axis, align.col.by, dot.size
) {
    if (is.null(align.col.by)) {
        colRange <- unique(dat[[col.axis]])
        colRange <- colRange[order(colRange)]
        colRange <- as.character(colRange)
    } else {
        colRange <- unique(align.col.by[[col.axis]])
        colRange <- colRange[order(colRange)]
        colRange <- as.character(colRange)
    }

    p <- ggplot2::ggplot(
        data = dat,
        aes(
            x = .data[[x.axis]], y = .data[[y.axis]], 
            colour = as.factor(.data[[col.axis]])
        )
    ) +
        ggplot2::geom_point(size = dot.size)
    if (!is.null(colRange)) {
        p <- p + ggplot2::lims(colour = colRange)
    }

    return(p)
}
#' Internal: Create a continuous colour scatter plot
#'
#' Generates a scatter or hexbin plot for continuous colouring of points.
#'
#' @param dat data.table containing the data to plot.
#' @param x.axis Character. Column name for x axis.
#' @param y.axis Character. Column name for y axis.
#' @param col.axis Character. Column name for colour (continuous).
#' @param colours Character. Colour scheme to use (e.g., "spectral", "jet", "viridis", etc).
#' @param dot.size Numeric. Size of dots.
#' @param hex Logical. Whether to use hex binning.
#' @param hex.bins Integer. Number of hex bins.
#' @param align.col.by a data.table. Used to align colour scale.
#' @param col.min.threshold Numeric. Minimum quantile for colour scale.
#' @param col.max.threshold Numeric. Maximum quantile for colour scale.
#'
#' @return ggplot2 object
#' @keywords internal
#' @noMd 
#' 
#' @import ggplot2
#' @import scales
#'
.make_continuous_scatter_plot <- function(dat, x.axis, y.axis, col.axis, colours, dot.size,
                                         hex, hex.bins, align.col.by, col.min.threshold, col.max.threshold) {

    
    colour.scheme <- .get_colour_scheme(colours)

    if (is.null(align.col.by)) {
        ColrMin <- quantile(dat[[col.axis]], probs = c(col.min.threshold), na.rm = TRUE)
        ColrMax <- quantile(dat[[col.axis]], probs = c(col.max.threshold), na.rm = TRUE)
    } else {
        ColrMin <- quantile(align.col.by[[col.axis]], probs = c(col.min.threshold), na.rm = TRUE)
        ColrMax <- quantile(align.col.by[[col.axis]], probs = c(col.max.threshold), na.rm = TRUE)
    }

    p <- ggplot2::ggplot(
        data = dat,
        aes(
            x = .data[[x.axis]],
            y = .data[[y.axis]],
            colour = .data[[col.axis]]
        )
    )

    if (hex) {
        p <- p + ggplot2::stat_summary_hex(
            aes(z = .data[[col.axis]]),
            fun = "mean",
            bins = hex.bins
        )
        p <- p + ggplot2::scale_fill_gradientn(
            colours = colour.scheme(50),
            limits = c(ColrMin, ColrMax),
            oob = scales::squish
        )
    } else {
        p <- p + ggplot2::geom_point(size = dot.size)
        p <- p + ggplot2::scale_colour_gradientn(
            colours = colour.scheme(50),
            limits = c(ColrMin, ColrMax),
            oob = scales::squish,
            na.value = "grey50"
        )
    }
    return(p)
}

#' Internal: Calculate centroids for factor groups using data.table
#'
#' @param dat data.table or data.frame
#' @param x.axis character, column name for x axis
#' @param y.axis character, column name for y axis
#' @param col.axis character, column name for grouping (factor)
#'
#' @return data.table with columns: centroidX, centroidY, centroidCol
#' @keywords internal
#' @import data.table
#' @noMd 
#' 
.calc_centroids_dt <- function(dat, x.axis, y.axis, col.axis) {
    centroids <- dat[, lapply(.SD, median), .SDcols = c(x.axis, y.axis), by = col.axis]
    setnames(centroids, x.axis, "centroidX")
    setnames(centroids, y.axis, "centroidY")
    setnames(centroids, col.axis, "centroidCol")
    return(centroids)
}

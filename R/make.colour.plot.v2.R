#' Make colour plot v2
#' 
#' Experimental feature.
#' Use at your own risk! 
#' 
#' This will replace the original colour plot in the future.
#' The idea is to refactor the original colour, fix up lingering bugs, and 
#' clean up the code in general.
#'
#' @author
#' Givanna Putri
#'
#'
#' @export

make.colour.plot <- function(dat,
                             x.axis,
                             y.axis,
                             col.axis = NULL,
                             col.type = "continuous",
                             # can be "continuous" or "factor"
                             add.label = FALSE,
                             # only works for 'factor'
                             hex = FALSE,
                             hex.bins = 30,
                             colours = "spectral",
                             # can be spectral, jet, only works for continuous 
                             col.min.threshold = 0.01,
                             col.max.threshold = 0.995,
                             align.xy.by = dat,
                             align.col.by = dat,
                             regression.line = NULL,
                             # "lm" # "loess"
                             title = col.axis,
                             filename = NULL,
                             dot.size = 1,
                             plot.width = 9,
                             plot.height = 7,
                             nudge_x = 0.5,
                             nudge_y = 0.5,
                             square = TRUE,
                             legend.loc = "right",
                             save.to.disk = TRUE,
                             path = getwd(),
                             blank.axis = FALSE,
                             col.tab = NULL) {
    ### Check for packages
    required_packages <- c(
        "ggplot2", "scales", "ggthemes", "viridis", "RColorBrewer"
    )
    check_packages_installed(required_packages)
    
    ### Load packages
    lapply(required_packages, function(x) require(x))
    
    
    ### Some tests
    
    if (hex == TRUE) {
        if (is.null(col.axis)) {
            message(
                "Note: hex bins do not currently work for density plots, 
                only for colour plots when col.axis is specified and can be 
                plotted as a continuous numeric variable"
            )
        }
        
        if (!is.null(col.axis)) {
            if (!is.numeric(dat[[col.axis]])) {
                # tests to see if there are any non numeric values
                stop(
                    "Sorry, hex bins only work when col.type is specified, 
                    and can be plotted as a continuous numeric variable"
                )
            }
        }
    }
    
    if (!is.null(col.axis)) {
        if (col.type == "continuous") {
            if (!is.numeric(dat[[col.axis]])) {
                # tests to see if there are any non numeric values
                message("Non-numeric values detected in col.axis 
                        -- using col.type = 'factor'")
                col.type <- "factor"
            }
        }
        
        if (col.type == "factor") {
            if (length(unique(as.factor(dat[[col.axis]]))) > 200) {
                message(
                    "Over 200 factors detected, using continuous scale 
                    instead of a factor scale"
                )
                col.type <- "continuous"
            }
        }
    }
    
    
    ### Define limits
    
    # X AXIS
    if (is.null(align.xy.by)) {
        Xmax <- max(dat[[x.axis]])
        Xmin <- min(dat[[x.axis]])
    } else {
        Xmax <- max(align.xy.by[[x.axis]])
        Xmin <- min(align.xy.by[[x.axis]])
    }
    
    # Y AXIS
    if (is.null(align.xy.by)) {
        Ymax <- max(dat[[y.axis]])
        Ymin <- min(dat[[y.axis]])
    } else {
        Ymax <- max(align.xy.by[[y.axis]])
        Ymin <- min(align.xy.by[[y.axis]])
    }
    
    # COLOUR
    
    if (!is.null(col.axis)) {
        if (col.type == "continuous") {
            if (is.null(align.col.by)) {
                ColrMin <-
                    quantile(dat[[col.axis]],
                             probs = c(col.min.threshold),
                             na.rm = TRUE)
                ColrMax <-
                    quantile(dat[[col.axis]],
                             probs = c(col.max.threshold),
                             na.rm = TRUE)
            } else {
                ColrMin <-
                    quantile(
                        align.col.by[[col.axis]],
                        probs = c(col.min.threshold),
                        na.rm = TRUE
                    )
                ColrMax <-
                    quantile(
                        align.col.by[[col.axis]],
                        probs = c(col.max.threshold),
                        na.rm = TRUE
                    )
            }
        }
        
        if (col.type == "factor") {
            if (is.null(align.col.by)) {
                
                colRange <- unique(dat[[col.axis]])
                colRange <- colRange[order(colRange)]
                colRange <- as.character(colRange)
            } else {
                
                colRange <- unique(align.col.by[[col.axis]])
                colRange <- colRange[order(colRange)]
                colRange <- as.character(colRange)
            }
        }
    }
    
    ### Initialise plot
    
    if (!is.null(col.axis)) {
        if (col.type == "continuous") {
            p <- ggplot(data = dat,
                        aes(
                            x = .data[[x.axis]],
                            y = .data[[y.axis]],
                            colour = .data[[col.axis]]
                        ))
            
            if (hex == TRUE) {
                p <- p + stat_summary_hex(aes(z = dat[[col.axis]]),
                                          fun = "mean",
                                          bins = hex.bins)
                p <- p + scale_fill_gradientn(
                    colours = c(colour.scheme(50)),
                    limits = c(ColrMin, ColrMax),
                    oob = squish
                )
            } else {
                p <- p + geom_point(size = dot.size)
                p <- p + scale_colour_gradientn(
                    colours = colour.scheme(50),
                    limits = c(ColrMin, ColrMax),
                    oob = squish,
                    na.value = "grey50"
                )
            }
        } else if (col.type == "factor") {
            if ("data.table" %in% class(colours)) {
                p <- ggplot(data = dat,
                            aes(
                                x = .data[[x.axis]],
                                y = .data[[y.axis]],
                                colour = as.factor(.data[[col.axis]])
                            )) +
                    geom_point(size = dot.size) +
                    lims(colour = colRange)
                
                col.vc <- setNames(as.character(colours[[1]]), colours[[2]])
                col.vc
                
                p <- p + scale_colour_manual(values = group.colors)
            } else {
                p <- ggplot(data = dat,
                            aes(
                                x = .data[[x.axis]],
                                y = .data[[y.axis]],
                                colour = as.factor(.data[[col.axis]])
                            )) +
                    geom_point(size = dot.size) +
                    lims(colour = colRange)
            }
        }
    }
    
    if (is.null(col.axis)) {
        p <- ggplot(data = dat,
                    aes(x = .data[[x.axis]],
                        y = .data[[y.axis]])) +
            ggpointdensity::geom_pointdensity(size = dot.size)
        
        if (colours == "viridis" ||
            colours == "magma" || colours == "inferno") {
            p <- p + viridis::scale_colour_viridis(option = colours)
        } else if (colours == "jet") {
            p <-
                p + ggplot2::scale_colour_gradientn(
                    colours = c(
                        "#00007F",
                        "blue",
                        "#007FFF",
                        "cyan",
                        "#7FFF7F",
                        "yellow",
                        "#FF7F00",
                        "red",
                        "#7F0000"
                    )
                )
        } else if (colours == "spectral") {
            p <-
                p + ggplot2::scale_colour_gradientn(colours = rev(
                    colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(50)
                ))
        }
        
        # Blue to Purple
        else if (colours == "BuPu") {
            colour.list <-
                (colorRampPalette(RColorBrewer::brewer.pal(9, "BuPu"))(31)) # 256
            # colours <- colorRampPalette(c(colour.list))
            p <-
                p + ggplot2::scale_colour_gradientn(colours = colour.list)
        }
    }
    
    ### Regression lione
    
    if (!is.null(regression.line)) {
        p <- p + geom_smooth(method = regression.line)
    }
    
    ### Add title
    
    if (is.null(title)) {
        title <- "Density"
    }
    
    p <- p + ggtitle(title)
    
    ### Set up axis
    p <-
        p + scale_x_continuous(
            breaks = scales::pretty_breaks(n = 8),
            name = x.axis,
            limits = c(Xmin, Xmax)
        )
    p <-
        p + scale_y_continuous(
            breaks = scales::pretty_breaks(n = 8),
            name = y.axis,
            limits = c(Ymin, Ymax)
        )
    
    ### Set up themes etc
    
    if (col.type == "continuous") {
        p <- p + theme(
            panel.background = element_rect(
                fill = "white",
                colour = "black",
                size = 0.5
            ),
            # change 'colour' to black for informative axis
            axis.title.x = element_text(color = "Black",
                                        # face="bold",
                                        size = 28),
            axis.title.y = element_text(color = "Black",
                                        # face="bold",
                                        size = 28),
            
            # axis.ticks = element_line(size = 4),
            # axis.line = element_line(size = 1),
            axis.text.x = element_text(color = "Black", size = 24),
            axis.text.y = element_text(color = "Black", size = 24),
            panel.border = element_rect(
                colour = "black",
                fill = NA,
                size = 2
            ),
            plot.title = element_text(
                color = "Black",
                face = "bold",
                size = 32,
                hjust = 0
            ) # size 70 for large, # 18 for small
        )
    }
    
    if (col.type == "factor") {
        p <- p + theme(
            panel.background = element_rect(
                fill = "white",
                colour = "black",
                size = 0.5
            ),
            axis.title.x = element_text(color = "Black",
                                        size = 28),
            axis.title.y = element_text(color = "Black",
                                        size = 28),
            
            axis.text.x = element_text(color = "Black", size = 24),
            axis.text.y = element_text(color = "Black", size = 24),
            panel.border = element_rect(
                colour = "black",
                fill = NA,
                size = 2
            ),
            plot.title = element_text(
                color = "Black",
                face = "bold",
                size = 32,
                hjust = 0
            ) # size 70 for large, # 18 for small
        )
        
        # p <- p + theme(legend.position="bottom")
    }
    
    if (square == TRUE) {
        p <- p + theme(aspect.ratio = 1)
    }
    
    ### Setup legend
    
    ## 'top' or 'bottom'
    if (legend.loc %in% c("top", "bottom")) {
        p <- p + theme(
            legend.direction = "horizontal",
            legend.position = legend.loc,
            legend.text = element_text(size = 18),
            legend.title = element_blank()
        )
    }
    
    ## 'left' or 'right'
    if (legend.loc %in% c("left", "right")) {
        p <- p + theme(
            legend.direction = "vertical",
            legend.position = legend.loc,
            legend.text = element_text(size = 18),
            legend.title = element_blank()
        )
    }
    
    ### Add labels (if desired)
    
    if (col.type == "factor") {
        if (add.label == TRUE) {
            ## Prepare centroids
            if (is.numeric(dat[[col.axis]])) {
                centroidX <- tapply(dat[[x.axis]], dat[[col.axis]], median) # median
                centroidY <- tapply(dat[[y.axis]], dat[[col.axis]], median)
                centroidCol <-
                    tapply(dat[[col.axis]], dat[[col.axis]], median)
                
                centroidsDf <- data.frame(centroidX, centroidY, centroidCol)
            }
            
            if (!is.numeric(dat[[col.axis]])) {
                labels <- sort(unique(dat[[col.axis]]))
                
                centroidsDf <- data.frame(
                    centroidX = tapply(dat[[x.axis]], dat[[col.axis]], median),
                    # median
                    centroidY = tapply(dat[[y.axis]], dat[[col.axis]], median),
                    centroidCol = labels
                )
            }
            
            ## Add labels
            p <- p + geom_point(
                data = centroidsDf,
                aes(x = centroidX,
                    y = centroidY),
                col = "black",
                # shape = 1,
                size = 2
            )
            
            p <- p + geom_label(
                data = centroidsDf,
                hjust = 0,
                nudge_x = nudge_x,
                nudge_y = nudge_y,
                aes(
                    x = centroidX,
                    y = centroidY,
                    label = centroidCol,
                    alpha = 0.5
                ),
                col = "black",
                fontface = "bold"
            )
            
            p <- p + guides(alpha = "none")
        }
    }
    
    ### Blank axis options
    
    if (blank.axis == TRUE) {
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
            # legend.position = "right",
            # legend.text=element_text(size=15), # large = 30 # small = 8
            # legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
            # legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
            # legend.title=element_blank(),
            # plot.title = element_text(color="Black", face="bold", size=15, hjust=0
        )
    }
    
    ### Save plot
    if (save.to.disk == TRUE) {
        if (!is.null(col.axis)) {
            if (col.type == "continuous") {
                lb <- "Colour"
            }
            
            if (col.type == "factor") {
                lb <- "Factor"
            }
        }
        
        if (is.null(col.axis)) {
            lb <- "Density plot"
        }
        
        if (is.null(filename)) {
            filename <-
                paste0(lb,
                       " plot - ",
                       title,
                       " - plotted on ",
                       x.axis,
                       " by ",
                       y.axis,
                       ".png")
        }
        
        ggsave(
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
    
    ### Print plot
    # print(p)
    # maybe return, i'm not sure.
    return(p)
}

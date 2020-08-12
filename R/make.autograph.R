#' make.autograph - Creates plots using ggplot2
#'
#' This function allows you to generate quantitative plots of your data.
#' Creates scatter plots.
#' If interested in generating other plots (bar graphs) or adding statistics, we recommend reading more on using 'ggplot2' (for more info https://ggplot2.tidyverse.org/).
#' Uses the package 'ggplot2' to generate plots, 'data.table' to manipulate data.
#'
#' @param dat NO DEFAULT. data.frame.
#' @param x.axis NO DEFAULT. Character or numeric. Selects column that will specify x-axis.
#' @param y.axis NO DEFAULT. Character or numeric. Selects column that will specify y-axis.
#' @param colour.by NO DEFAULT. Character or numeric. Selects column that will determine differentiator of colour in plots.
#' @param colours NO DEFAULT. String of characters. Will determine the colours that will differentiate between colours specified by 'colour.by'.
#' @param y.axis.label NO DEFAULT. Character. Specify y-axis label. Can be "" if one is not desired.
#' @param grp.order DEFAULT = NULL. String of characters or numbers. Specify order of group for x-axis.
#' @param title DEFAULT = paste0(y.axis.label, " - ", y.axis). Character. Specify title of plot on graph. Default combines 'y.axis.label' and 'y.axis'.
#' @param filename DEFAULT = paste0(y.axis.label, " - ", y.axis, ".pdf"). Character. Specify name plot will be exported as. Default combines 'y.axis.label' and 'y.axis'. Make sure to end with ".pdf" so file is correctly saved.
#' @param scale DEFAULT = "lin". Character. Select y-axis scale. Can also be "sci".
#' @param dot.size DEFAULT = 5. Numeric. Specify size of dots on plot.
#' @param width DEFAULT = 5. Numeric. Specify width of plot.
#' @param height DEFAULT = 5. Numeric. Specify height of plot.
#' @param path DEFAULT = getwd(). The location to save plots. By default, will save to current working directory. Can be overidden.
#'
#' @usage make.autograph(dat, x.axis, y.axis, colour.by, colours, y.axis.label, grp.order, title, filename, scale, dot.size, width, height, path)
#'
#' @examples
#' dat <- data.frame(Samples = c("Mock_01", "Mock_02", "Mock_03", "WNV_01", "WNV_02", "WNV_03"),
#'                   Group = c(rep("Mock", 3), rep("WNV", 3)),
#'                   Tcells = c(20, 40, 30, 60, 70, 80),
#'                   Bcells = c(90, 95, 70, 20, 15, 30),
#'                   Batch = c(1,2,1,3,2,1)
#'                   )
#'
#' Spectre::make.autograph(dat = dat,
#'                         x.axis = "Group",
#'                         y.axis = "Tcells",
#'                         colour.by = "Batch",
#'                         colours = c("Black", "Red", "Blue"),
#'                         y.axis.label = "Proportion"
#'                         )
#'
#' @author
#' Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#' @export

make.autograph <- function(dat,
                           x.axis,
                           y.axis,

                           colour.by = x.axis,

                           y.axis.label = y.axis,

                           grp.order = NULL,
                           colours = NULL,

                           my_comparisons = NULL,
                           Variance_test = NULL,
                           Pairwise_test = NULL,

                           title = paste0(y.axis),
                           subtitle = NULL,
                           filename = paste0(y.axis, ".pdf"),

                           violin = TRUE,
                           scale = "lin", # can be "lin" or "sci"
                           dot.size = 5,
                           width = 5,
                           height = 5,
                           max.y = 1.4,

                           path = getwd()
                           # my_comparisons = NULL,
                           # y.first.level = 1.2,
                           # y.second.level = 1.5,
                           # Variance_test = "kruskal.test", # or "anova", or NULL
                           # Pairwise_test = "wilcox.test", # or 't.test', or NULL
                           )
{
  ### Test data
      # if(!require('ggplot2')) {install.packages('ggplot2')}
      # if(!require('ggpubr')) {install.packages('ggpubr')}
      # if(!require('scales')) {install.packages('scales')}
      # if(!require('devtools')) {install.packages('devtools')}
      # if(!require('rstudioapi')) {install.packages('rstudioapi')}
      #
      # scale <- "lin"
      # dot.size <- 5
      # width <- 5
      # height <- 5
      # y.first.level = 1.2
      # y.second.level = 1.5
      # Variance_test <- "kruskal.test"
      # Pairwise_test <- "wilcox.test"

  ### Check that necessary packages are installed
      if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
      if(!is.element('ggplot2', installed.packages()[,1])) stop('ggplot2 is required but not installed')
      if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')

  ### Require packages
      require(Spectre)
      require(ggplot2)
      require(data.table)

  ### Library test
  if(!is.null(colours)){
    if(length(unique(dat[[colour.by]])) != length(colours)){
      stop('The length of factors you want to colour the plot by does not match the number of colours you have provided.')
    }
  }

  ### Set up colours and other settings

  message(paste0("AutoGraph for `", y.axis.label, " - ", y.axis, "` started"))

  message("AutoGraph - setup started")

    # spectral.list <- colorRampPalette(brewer.pal(11,"Spectral"))(50)
    # spectral.list <- rev(spectral.list)
    # colour.scheme <- colorRampPalette(c(spectral.list))

      spectral.list <- colorRampPalette(brewer.pal(11,"Spectral"))(50)
      spectral.list <- rev(spectral.list)
      colour.scheme <- colorRampPalette(c(spectral.list))

  ### Max and min values
    dat <- data.table::as.data.table(dat)

    max_y_value <- max(dat[,y.axis, with = FALSE], na.rm = TRUE)
    # max_y_value_p10 <- max_y_value*y.first.level
    # max_y_value_p40 <- max_y_value*y.second.level
    max_y_value_p40 <- max_y_value*max.y

    min_y_value<- min(dat[,y.axis, with = FALSE], na.rm = TRUE)
    bottom_y <- min_y_value
    bottom_y

  ###
    Xaxis <- dat[[x.axis]] <- as.factor(dat[[x.axis]])

    if(is.null(grp.order)){
      if(is.null(colours)){
        colours <- colour.scheme(length(unique(dat[[x.axis]])))
      }
    }

    if(!is.null(grp.order)){
      Xaxis <- dat[[x.axis]] <- factor(dat[[x.axis]], levels = grp.order)

      if(is.null(colours)){
        colours <- colour.scheme(length(as.vector(grp.order)))
      }
    }



  ###

  message("AutoGraph - setup complete")

  ### Plotting

  message("AutoGraph - plotting started")

      p <- ggplot2::ggplot(dat, ggplot2::aes(x=Xaxis, y=dat[[y.axis]])) #, fill=Species)) + # can use fill = Species -- makes coloured by Species, with block outline of point -- NEED FILL for violin plot to be filled

      ## POINTS
      p <- p + ggplot2::geom_point(ggplot2::aes(fill=as.factor(dat[[colour.by]]),
                                                colour = as.factor(dat[[colour.by]])),
                          shape = 21,
                          stroke = 0,
                          size = dot.size,
                          position = ggplot2::position_jitter(width = 0.1, height = 0))
      ## VIOLIN
      if(violin == TRUE){
        message("AutoGraph - adding violin plot")
        p <- p + geom_violin(ggplot2::aes(fill=as.factor(dat[[colour.by]]),
                                          colour = as.factor(dat[[colour.by]])),
                             trim=FALSE,
                             show.legend = FALSE,
                             alpha = 0.1)
      }


      ## COLOUR CONTROL
          p <- p + ggplot2::scale_fill_manual(name = colour.by, values = colours)
          p <- p + ggplot2::scale_color_manual(name = colour.by, values = colours)

          p <- p + ggplot2::ggtitle(title)
          p <- p + ggplot2::labs(x= paste0(x.axis),
                                 y = y.axis.label,
                                 subtitle = subtitle) # colnames(data)[3] would return the name -- use similar for loop -- maybe data$dose

      # OTHER OPTIONS
          #scale_fill_brewer(palette="Dark2") +
          #scale_color_brewer(palette="Dark2") +

      ## PRISM -- SE (from https://stackoverflow.com/questions/2676554/in-r-how-to-find-the-standard-error-of-the-mean)

            p <- p + ggplot2::stat_summary(
              fun.max=function(i) mean(i) + sd(i)/sqrt(length(i)), # 0.975
              fun.min=function(i) mean(i) - sd(i)/sqrt(length(i)), # 0.975
              geom="errorbar", width=0.5, size = 1)

            p <- p + ggplot2::stat_summary(fun=mean,
                                           fun.min = mean,
                                           fun.max = mean, geom="crossbar",
                                           width = 0.7,
                                           size = 0.5) # add a large point for the mean

      # VARIANCE
          #stat_summary(fun.y=mean, geom="point", shape=18, size=8, color="Black") + # add a large point for the mean
          #geom_pointrange() + # need ymin and ymax
          #geom_errorbar(aes(ymax = Sepal.Length+se, ymax = Sepal.Length-se)) + # SEM -- requires prior calculation of SE
          #geom_errorbar(aes(ymax = Sepal.Length+sd, ymax = Sepal.Length-sd)) +

      # MORE THAN TWO GROUPS: pairwise comparison with overall anova/Kruskal-Wallis result

          if(!is.null(my_comparisons)){
            if(!is.null(Pairwise_test)){
              p <- p + stat_compare_means(comparisons = my_comparisons, method = Pairwise_test) #, label.y = max_y_value_p10) + # Add pairwise comparisons p-value # default is "wilcox.test" (non-parametric), can be "t.test" (parametric) # Add global Anova ("anova") or Kruskal-Wallis ("kruskal.test", default) p-value # an add label.y = 50 to specifiy position
            }
          }

          if(!is.null(Variance_test)){
            p <- p + stat_compare_means(method = Variance_test, label.y = max_y_value_p40, size = 4) # Add global Anova ("anova") or Kruskal-Wallis ("kruskal.test", default) p-value # an add label.y = 50 to specifiy position
          }

      # MORE THAN TWO GROUPS: compare against reference sample
          #stat_compare_means(method = "kruskal.test", label.y = 45) +      # Add global p-value
          #stat_compare_means(label = "p.signif", method = "t.test",ref.group = "0.5", label.y = 40)

      ## VISUALS
      #labs(title=paste0(i, " ", plotname), x= paste0(X_axis_label), y = paste0(Y_axis_label)) + # colnames(data)[3] would return the name -- use similar for loop -- maybe data$dose


      #coord_fixed(ratio = 1) + # determines size ratio of the plot -- smaller increases width
          p <- p + ggplot2::theme_classic(base_size = 30) # can be theme_classic(), theme_bw()

      ## THEMES
          p <- p + ggplot2::theme(legend.position = "right", # can be "left" "right" "top" "bottom" "none
                legend.text = ggplot2::element_text(colour="black",size=10,angle=0,hjust=0,vjust=0,face="bold"),
                legend.title = ggplot2::element_text(colour="black",size=10,angle=0,hjust=0,vjust=0,face="bold"),

                axis.text.x = ggplot2::element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="bold"),
                axis.text.y = ggplot2::element_text(colour="black",size=12,angle=0,hjust=1,vjust=0,face="bold"),
                axis.title.x = ggplot2::element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0,face="bold"),
                axis.title.y = ggplot2::element_text(colour="black",size=12,angle=90,hjust=.5,vjust=1,face="bold"),

                # Title
                plot.title = ggplot2::element_text(lineheight=.8, face="bold", hjust = 0, size = 12), # hjust = 0.5 to centre
                plot.subtitle = ggplot2::element_text(size=10, #åhjust=0.5,
                                                      face="italic", color="black"),
                axis.line = ggplot2::element_line(colour = 'black', size = 1),
                axis.ticks = ggplot2::element_line(colour = "black", size = 1)
                )

    ### End construction of 'p'

    ## SCALES
    #scale_x_discrete(limits=Group_Order) + # use to re-arrange X axis values
    if(scale == "lin"){
      p <- p + ggplot2::scale_y_continuous(limits = c(0, max_y_value_p40))
    }

    if(scale == "sci"){
      p <- p + ggplot2::scale_y_continuous(labels = scales::scientific, limits = c(0, max_y_value_p40))
    }

    # if(scale == "log2"){
    #   p <- p + scale_y_continuous(trans = "log2")
    # }
    #
    # if(scale == "log10"){
    #   p <- p + scale_y_continuous(trans = "log10")
    #   # p <- p + scale_y_log10()
    # }

    ## View plot
    p

    ## Save
    ggplot2::ggsave(plot = p, filename = paste0(filename), width = width, height = height, path = path) # width 3.6 default, height 5 default
    message(paste0("AutoGraph for `", y.axis.label, " - ", y.axis, "` saved to disk"))
}

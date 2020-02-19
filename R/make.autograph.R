#' make.autograph
#'
#' @export

make.autograph <- function(x,
                           x.axis,
                           y.axis,
                           colour.by,

                           colours,
                           y.axis.label,

                           my_comparisons,

                           title = paste0(y.axis.label, " - ", y.axis),
                           filename = paste0(y.axis.label, " - ", y.axis, ".pdf"),
                           path = getwd(),

                           scale = "lin",
                           dot.size = 5,
                           width = 5,
                           height = 5,
                           y.first.level = 1.2,
                           y.second.level = 1.5,
                           Variance_test = "kruskal.test",
                           Pairwise_test = "wilcox.test"
                           )
{
  ### Test data
      # if(!require('ggplot2')) {install.packages('ggplot2')}
      # if(!require('ggpubr')) {install.packages('ggpubr')}
      # if(!require('scales')) {install.packages('scales')}
      # if(!require('devtools')) {install.packages('devtools')}
      # if(!require('rstudioapi')) {install.packages('rstudioapi')}
      #
      # library("Spectre")
      # library(ggplot2)
      # library(ggpubr)
      # library(scales)
      # library(devtools)
      # library(rstudioapi)
      #
      # x <- data.frame(Samples = c("Mock_01", "Mock_02", "Mock_03", "WNV_01", "WNV_02", "WNV_03"),
      #                 Groups = c(rep("Mock", 3), rep("WNV", 3)),
      #                 Tcells = c(200, 400, 300, 600, 700, 800),
      #                 Bcells = c(900, 950, 700, 200, 150, 300),
      #                 Batches = c(1,2,1,3,2,1)
      #                 )
      #
      # x.axis <- "Groups"
      # y.axis <- "Tcells"
      # colour.by <- "Batches"
      #
      # colours <- c("Black", "Red", "Blue")
      # y.axis.label <- "Frequency"
      #
      # my_comparisons <- list(c("Mock", "WNV"))
      #
      # title <- paste0(y.axis.label, " - ", y.axis)
      # filename <-  paste0(y.axis.label, " - ", y.axis, ".pdf")
      #
      # scale <- "lin"
      # dot.size <- 5
      # width <- 5
      # height <- 5
      # y.first.level = 1.2
      # y.second.level = 1.5
      # Variance_test <- "kruskal.test"
      # Pairwise_test <- "wilcox.test"

  ### Library test
      if(length(unique(x[[colour.by]])) != length(colours)) stop('The length of factors you want to colour the plot by does not match the number of colours you have provided.')

  ### Set up colours and other settings
    # spectral.list <- colorRampPalette(brewer.pal(11,"Spectral"))(50)
    # spectral.list <- rev(spectral.list)
    # colour.scheme <- colorRampPalette(c(spectral.list))

  ### Max and min values
    max_y_value <- max(x[y.axis], na.rm = TRUE)
    max_y_value_p10 <- max_y_value*y.first.level
    max_y_value_p40 <- max_y_value*y.second.level

    min_y_value<- min(x[y.axis], na.rm = TRUE)
    bottom_y <- min_y_value
    bottom_y

  ###
    Xaxis <- x[[x.axis]] <- as.factor(x[[x.axis]])

  ### Plotting

    p <- ggplot(x, aes(x=Xaxis, y=x[[y.axis]])) #, fill=Species)) + # can use fill = Species -- makes coloured by Species, with block outline of point -- NEED FILL for violin plot to be filled

      ## POINTS
          p <- p + geom_point(aes(fill=as.factor(x[[colour.by]]), colour = as.factor(x[[colour.by]])), shape=21, stroke = 0, size = dot.size, position=position_jitter(width = 0.1, height = 0))

      ## COLOUR CONTROL
          p <- p + scale_fill_manual(name = colour.by, values = colours)
          p <- p + scale_color_manual(name = colour.by, values = colours)

          p <- p + ggtitle(title)
          p <- p + labs(x= paste0(x.axis), y = y.axis.label) # colnames(data)[3] would return the name -- use similar for loop -- maybe data$dose

      # OTHER OPTIONS
          #scale_fill_brewer(palette="Dark2") +
          #scale_color_brewer(palette="Dark2") +

      ## PRISM -- SE (from https://stackoverflow.com/questions/2676554/in-r-how-to-find-the-standard-error-of-the-mean)
          p <- p + stat_summary(
                  fun.ymax=function(i) mean(i) + sd(i)/sqrt(length(i)), # 0.975
                  fun.ymin=function(i) mean(i) - sd(i)/sqrt(length(i)), # 0.975
                  geom="errorbar", width=0.5, size = 1)

          p <- p + stat_summary(fun.y=mean, fun.ymin = mean, fun.ymax = mean, geom="crossbar", width = 0.7, size = 0.5) # add a large point for the mean


      # VARIANCE
          #stat_summary(fun.y=mean, geom="point", shape=18, size=8, color="Black") + # add a large point for the mean
          #geom_pointrange() + # need ymin and ymax
          #geom_errorbar(aes(ymax = Sepal.Length+se, ymax = Sepal.Length-se)) + # SEM -- requires prior calculation of SE
          #geom_errorbar(aes(ymax = Sepal.Length+sd, ymax = Sepal.Length-sd)) +

      # MORE THAN TWO GROUPS: pairwise comparison with overall anova/Kruskal-Wallis result
          p <- p + stat_compare_means(comparisons = my_comparisons, method = Pairwise_test) #, label.y = max_y_value_p10) + # Add pairwise comparisons p-value # default is "wilcox.test" (non-parametric), can be "t.test" (parametric)
          p <- p + stat_compare_means(method = Variance_test, label.y = max_y_value_p40, size = 4) # Add global Anova ("anova") or Kruskal-Wallis ("kruskal.test", default) p-value # an add label.y = 50 to specifiy position

      # MORE THAN TWO GROUPS: compare against reference sample
          #stat_compare_means(method = "kruskal.test", label.y = 45) +      # Add global p-value
          #stat_compare_means(label = "p.signif", method = "t.test",ref.group = "0.5", label.y = 40)

      ## VISUALS
      #labs(title=paste0(i, " ", plotname), x= paste0(X_axis_label), y = paste0(Y_axis_label)) + # colnames(data)[3] would return the name -- use similar for loop -- maybe data$dose


      #coord_fixed(ratio = 1) + # determines size ratio of the plot -- smaller increases width
          p <- p + theme_classic(base_size = 30) # can be theme_classic(), theme_bw()

      ## THEMES
          p <- p + theme(legend.position = "right", # can be "left" "right" "top" "bottom" "none
                legend.text = element_text(colour="black",size=10,angle=0,hjust=0,vjust=0,face="bold"),
                legend.title = element_text(colour="black",size=10,angle=0,hjust=0,vjust=0,face="bold"),

                axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="bold"),
                axis.text.y = element_text(colour="black",size=12,angle=0,hjust=1,vjust=0,face="bold"),
                axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0,face="bold"),
                axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=1,face="bold"),

                # Title
                plot.title = element_text(lineheight=.8, face="bold", hjust = 0, size = 12), # hjust = 0.5 to centre
                axis.line = element_line(colour = 'black', size = 1),
                axis.ticks = element_line(colour = "black", size = 1)
                )

    ### End construction of 'p'

    ## SCALES
    #scale_x_discrete(limits=Group_Order) + # use to re-arrange X axis values
    if(scale == "lin"){
      p <- p + scale_y_continuous(limits = c(0, max_y_value_p40))
      }
    if(scale == "sci"){
      p <- p + scale_y_continuous(labels = scales::scientific, limits = c(0, max_y_value_p40))
      }

    ## View plot
    p

    ## Save
    ggsave(plot = p, filename = paste0(filename), width = 5, height = 5, path = path) # wdith 3.6 default, height 5 default

}

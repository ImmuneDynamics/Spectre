# Spectre::colour.plot

colour.plot <- function(d,
                        x.axis, # "UMAP1"
                        y.axis, # "UMAP2"
                        col.axis, # "BV605.Ly6C"
                        col.min.threshold = 0.01,
                        col.max.threshold = 0.995,
                        title,
                        colours,
                        dot.size,
                        align.xy.by = NULL, # choose a data frame to set absolute limits for X/Y/colour
                        align.col.by = NULL
                        # ### align axis and colours
                            # align.x = FALSE,
                            # align.y = FALSE,
                            # align.colour = FALSE,

                        ## Manual set max and min for alignment
                            # xMax = NULL,
                            # xMin = NULL,
                            #
                            # yMax = NULL,
                            # yMin = NULL,
                            #
                            # colourMax = NULL,
                            # colourMin = NULL
                        ){

  ## -- copy the # codes for 'spectral' colours -- use that directly, avoid needing Rcolourbrewer etc

  ## TESTING
      # d <- plot.dat
      # x.axis = "UMAP_42_X"
      # y.axis = "UMAP_42_Y"
      # col.axis = "BV605.Ly6C"
      # col.min.threshold = 0.01
      # col.max.threshold = 1.0
      # title = paste0("All samples", " - ", "BV605.Ly6C")
      # colours = "spectral"
      # dot.size = 1
      #
      # align.xy.by = plot.dat
      # align.col.by = plot.dat
      #
      # xMax = 5
      # xMin = -5
      # yMax = 5
      # yMin = -5
      # colourMax = 600
      # colourMin = 0


  ## Colour setup
      if(colours == "jet"){
        jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
        colour.scheme <- jet
      }

      if(colours == "spectral"){
        spectral.list <- colorRampPalette(brewer.pal(11,"Spectral"))(50)
        spectral.list <- rev(spectral.list)
        spectral <- colorRampPalette(c(spectral.list))
        colour.scheme <- spectral
      }


  ## Define limits

      # X AXIS
      if(is.null(align.xy.by) == TRUE){
        Xmax <- max(d[[x.axis]])
        Xmin <- min(d[[x.axis]])
      }

      if(is.null(align.xy.by) == FALSE){
        Xmax <- max(align.xy.by[[x.axis]])
        Xmin <- min(align.xy.by[[x.axis]])
        }


      # Y AXIS
      if(is.null(align.xy.by) == TRUE){
        Ymax <- max(d[[y.axis]])
        Ymin <- min(d[[y.axis]])
        }

      if(is.null(align.xy.by) == FALSE){
        Ymax <- max(align.xy.by[[y.axis]])
        Ymin <- min(align.xy.by[[y.axis]])
        }

      # COLOUR
      if(is.null(align.col.by) == TRUE){
        ColrMin <- quantile(d[[col.axis]], probs = c(col.min.threshold))
        ColrMax <- quantile(d[[col.axis]], probs = c(col.max.threshold))
      }

      if(is.null(align.col.by) == FALSE){
        ColrMin <- quantile(align.col.by[[col.axis]], probs = c(col.min.threshold))
        ColrMax <- quantile(align.col.by[[col.axis]], probs = c(col.max.threshold))
      }


  ## Generate and show coloured plot
  print(ggplot(data = d, aes(x = d[[x.axis]], y = d[[y.axis]], colour = d[[col.axis]])) +
    geom_point(size = dot.size) +

    #scale_colour_gradientn(colours = c("Black", colour.scheme(50)),
    #                       #values = rescale(c(ColrMin, ColrMax)),
    #                       breaks = c(ColrMin),
    #                       limits = c(min(d[[col.axis]]), ColrMax),
    #                       na.value = "Purple"
    #                       )+

    scale_colour_gradientn(colours = colour.scheme(50),
                           limits = c(ColrMin, ColrMax), #0.03-01 seems to work well#0.97-995 seems to work well
                           #na.value = "Purple")+
                           oob=squish) +
    ggtitle(title) +
    #guides(fill = guide_legend(angle = 90, colour = "red"))+ ###
    #labs(colour = col.axis)+

    xlim(Xmin, Xmax) +
    ylim(Ymin, Ymax) +

    xlab(x.axis)+
    ylab(y.axis)+
      theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # change 'colour' to black for informative axis
      #axis.line=element_blank(),
      #axis.text.x=element_blank(),
      #axis.text.y=element_blank(),
      #axis.ticks=element_blank(),
      axis.title.x=element_text(color="Black", face="bold", size=18),
      axis.title.y=element_text(color="Black", face="bold", size=18),
      #panel.grid.major = element_blank(),
      #panel.background=element_blank(),
      #panel.border=element_blank(),
      #panel.grid.minor=element_blank(),
      #plot.background=element_blank(),
      #legend.position = "right",
      legend.text=element_text(size=12), # large = 30 # small = 8
      legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
      legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
      legend.title=element_blank(),
      plot.title = element_text(color="Black", face="bold", size=22, hjust=0) # size 70 for large, # 18 for small
      )
  )

}
